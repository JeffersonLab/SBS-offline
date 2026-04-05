#!/usr/bin/env python3
"""
Generate AIML training data from G4SBS simulation ROOT files.

Reads G4SBS Monte Carlo truth trees and produces preprocessed training
datasets for the autoencoder denoiser and GNN track finder.

Supports the standard G4SBS output formats used at JLab:
- Raw G4SBS output (truth hits + tracks from Geant4)
- Digitized output (strip-level ADC waveforms after SBS-offline digitization)

The generator handles:
- Extraction of truth-level GEM hits with full MC particle information
- Labeling of signal vs background hits using track traversal criteria
- Construction of graph representations for GNN training
- Strip-level waveform extraction for convolutional autoencoder training
- Train/validation/test splitting with event-level separation
- Normalization statistics computation

Usage:
    # From G4SBS truth-level output (e.g., GMn simulation):
    python generate_training_data.py \\
        --input /volatile/halla/sbs/simulations/gmn_sbs4_*.root \\
        --det-prefix Earm_BBGEM \\
        --output-dir training_data/gmn \\
        --max-events 50000

    # From digitized output (includes ADC waveforms):
    python generate_training_data.py \\
        --input /path/to/digitized/simdigtest_*.root \\
        --det-prefix Earm_BBGEM \\
        --output-dir training_data/gmn_dig \\
        --include-waveforms \\
        --max-events 50000

    # Multiple files via glob:
    python generate_training_data.py \\
        --input '/volatile/halla/sbs/digitized/*.root' \\
        --output-dir training_data/combined

    # Quick test with synthetic data:
    python generate_training_data.py \\
        --synthetic --n-events 5000 \\
        --output-dir training_data/synthetic
"""

import argparse
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.root_io import (
    load_g4sbs_gem_data,
    load_g4sbs_digitized_strips,
    generate_synthetic_gem_data,
    hits_to_graph,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate AIML training data from G4SBS ROOT files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # GMn BigBite GEM truth data:
  python generate_training_data.py --input 'gmn_sbs4_*.root' --det-prefix Earm_BBGEM

  # With digitized strip waveforms:
  python generate_training_data.py --input 'simdigtest_*.root' --include-waveforms

  # Synthetic data for testing the pipeline:
  python generate_training_data.py --synthetic --n-events 5000
        """,
    )

    # Input sources
    input_group = parser.add_argument_group("Input")
    input_group.add_argument(
        "--input", type=str, nargs="+", default=None,
        help="Input ROOT file(s). Supports glob patterns.",
    )
    input_group.add_argument(
        "--synthetic", action="store_true",
        help="Generate synthetic training data instead of reading ROOT files.",
    )
    input_group.add_argument(
        "--n-events", type=int, default=10000,
        help="Number of synthetic events to generate (default: 10000).",
    )
    input_group.add_argument(
        "--tree-name", type=str, default="T",
        help="ROOT TTree name (default: T).",
    )
    input_group.add_argument(
        "--det-prefix", type=str, default="Earm_BBGEM",
        help="G4SBS detector branch prefix (default: Earm_BBGEM for BigBite).",
    )
    input_group.add_argument(
        "--max-events", type=int, default=None,
        help="Maximum events to read per file.",
    )

    # Output
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output-dir", type=str, default="training_data",
        help="Output directory for training data (default: training_data).",
    )
    output_group.add_argument(
        "--format", type=str, choices=["npz", "pt"], default="npz",
        help="Output format: npz (numpy) or pt (PyTorch). Default: npz.",
    )

    # Processing options
    proc_group = parser.add_argument_group("Processing")
    proc_group.add_argument(
        "--include-waveforms", action="store_true",
        help="Also extract digitized strip ADC waveforms.",
    )
    proc_group.add_argument(
        "--build-graphs", action="store_true", default=True,
        help="Pre-build graph representations for GNN training (default: True).",
    )
    proc_group.add_argument(
        "--no-build-graphs", action="store_false", dest="build_graphs",
        help="Skip pre-building graphs.",
    )
    proc_group.add_argument(
        "--r-connect", type=float, default=50.0,
        help="Edge construction radius in mm (default: 50).",
    )
    proc_group.add_argument(
        "--k-nearest", type=int, default=8,
        help="Max nearest neighbors per hit for edges (default: 8).",
    )
    proc_group.add_argument(
        "--min-hits-per-track", type=int, default=3,
        help="Minimum hits for a particle to be labeled as signal (default: 3).",
    )
    proc_group.add_argument(
        "--primary-only", action="store_true",
        help="Only label primary particles (trid==1) as signal.",
    )

    # Splitting
    split_group = parser.add_argument_group("Train/Val/Test Split")
    split_group.add_argument(
        "--train-fraction", type=float, default=0.7,
        help="Training set fraction (default: 0.7).",
    )
    split_group.add_argument(
        "--val-fraction", type=float, default=0.15,
        help="Validation set fraction (default: 0.15).",
    )
    split_group.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for splitting (default: 42).",
    )

    # Synthetic data options
    synth_group = parser.add_argument_group("Synthetic Data Options")
    synth_group.add_argument(
        "--n-layers", type=int, default=5,
        help="Number of GEM layers for synthetic data (default: 5).",
    )
    synth_group.add_argument(
        "--n-signal-tracks", type=int, default=2,
        help="Signal tracks per event for synthetic data (default: 2).",
    )
    synth_group.add_argument(
        "--n-noise-hits", type=int, default=10,
        help="Noise hits per layer for synthetic data (default: 10).",
    )

    return parser.parse_args()


def expand_input_files(input_patterns):
    """Expand glob patterns in input file list."""
    files = []
    for pattern in input_patterns:
        expanded = sorted(glob.glob(pattern))
        if not expanded:
            print(f"  Warning: no files match pattern '{pattern}'")
        files.extend(expanded)
    return files


def relabel_events(events, min_hits_per_track=3, primary_only=False):
    """Re-apply signal/noise labeling with custom criteria.

    Parameters
    ----------
    events : list[dict]
        Events loaded from G4SBS.
    min_hits_per_track : int
        Minimum hits across planes for a track to be labeled signal.
    primary_only : bool
        If True, only label primary particle (trid==1) as signal.
    """
    for event in events:
        trid = event["track_ids"]
        labels = np.zeros(len(trid), dtype=np.float32)

        if primary_only:
            labels[trid == 1] = 1.0
        else:
            unique_trids, counts = np.unique(trid[trid > 0], return_counts=True)
            signal_trids = set(unique_trids[counts >= min_hits_per_track])
            for i, tid in enumerate(trid):
                if int(tid) in signal_trids:
                    labels[i] = 1.0

        event["labels"] = labels

    return events


def compute_statistics(events):
    """Compute dataset-wide statistics for normalization and reporting."""
    all_hits = np.concatenate([ev["hits"] for ev in events if len(ev["hits"]) > 0])
    all_labels = np.concatenate([ev["labels"] for ev in events if len(ev["labels"]) > 0])

    n_total_hits = len(all_labels)
    n_signal = int(all_labels.sum())
    n_noise = n_total_hits - n_signal

    stats = {
        "n_events": len(events),
        "n_total_hits": n_total_hits,
        "n_signal_hits": n_signal,
        "n_noise_hits": n_noise,
        "signal_fraction": n_signal / n_total_hits if n_total_hits > 0 else 0,
        "hits_per_event_mean": n_total_hits / len(events) if events else 0,
        "hits_per_event_std": float(np.std([len(ev["hits"]) for ev in events])),
        "feature_means": all_hits.mean(axis=0).tolist(),
        "feature_stds": all_hits.std(axis=0).tolist(),
        "feature_mins": all_hits.min(axis=0).tolist(),
        "feature_maxs": all_hits.max(axis=0).tolist(),
    }

    # Track statistics
    n_tracks = []
    for ev in events:
        unique_trids = np.unique(ev["track_ids"][ev["track_ids"] > 0])
        n_tracks.append(len(unique_trids))
    stats["tracks_per_event_mean"] = float(np.mean(n_tracks)) if n_tracks else 0
    stats["tracks_per_event_std"] = float(np.std(n_tracks)) if n_tracks else 0

    # Layer statistics
    all_layers = np.concatenate([ev["layer_ids"] for ev in events])
    unique_layers = np.unique(all_layers)
    stats["n_layers"] = len(unique_layers)
    stats["layer_ids"] = unique_layers.tolist()
    stats["hits_per_layer"] = {
        int(l): int((all_layers == l).sum()) for l in unique_layers
    }

    return stats, all_hits.mean(axis=0), all_hits.std(axis=0)


def split_events(events, train_frac, val_frac, seed=42):
    """Split events into train/val/test sets."""
    rng = np.random.RandomState(seed)
    n = len(events)
    indices = rng.permutation(n)

    n_train = int(n * train_frac)
    n_val = int(n * val_frac)

    train_idx = indices[:n_train]
    val_idx = indices[n_train:n_train + n_val]
    test_idx = indices[n_train + n_val:]

    train = [events[i] for i in train_idx]
    val = [events[i] for i in val_idx]
    test = [events[i] for i in test_idx]

    return train, val, test


def save_split(events, output_path, fmt="npz", build_graphs=False,
               r_connect=50.0, k_nearest=8):
    """Save a set of events to disk."""
    if fmt == "npz":
        # Pack events into arrays
        all_hits = []
        all_layers = []
        all_labels = []
        all_track_ids = []
        event_boundaries = [0]

        for ev in events:
            all_hits.append(ev["hits"])
            all_layers.append(ev["layer_ids"])
            all_labels.append(ev["labels"])
            all_track_ids.append(ev["track_ids"])
            event_boundaries.append(event_boundaries[-1] + len(ev["hits"]))

        save_dict = {
            "hits": np.concatenate(all_hits) if all_hits else np.zeros((0, 14), dtype=np.float32),
            "layer_ids": np.concatenate(all_layers) if all_layers else np.zeros(0, dtype=np.int64),
            "labels": np.concatenate(all_labels) if all_labels else np.zeros(0, dtype=np.float32),
            "track_ids": np.concatenate(all_track_ids) if all_track_ids else np.zeros(0, dtype=np.int64),
            "event_boundaries": np.array(event_boundaries, dtype=np.int64),
        }

        # Optionally include pre-built graphs
        if build_graphs:
            all_edge_index = []
            all_edge_attr = []
            all_edge_labels = []
            graph_boundaries = [0]
            node_offset = 0

            for ev in events:
                if len(ev["hits"]) == 0:
                    graph_boundaries.append(graph_boundaries[-1])
                    continue

                edge_index, edge_attr = hits_to_graph(
                    ev["hits"], ev["layer_ids"],
                    r_connect=r_connect, k_nearest=k_nearest,
                )

                # Compute edge truth labels
                n_edges = edge_index.shape[1]
                edge_labels = np.zeros(n_edges, dtype=np.float32)
                for e in range(n_edges):
                    src, dst = edge_index[0, e], edge_index[1, e]
                    if (ev["labels"][src] == 1 and ev["labels"][dst] == 1 and
                            ev["track_ids"][src] > 0 and
                            ev["track_ids"][src] == ev["track_ids"][dst]):
                        edge_labels[e] = 1.0

                # Offset node indices for global graph
                edge_index_global = edge_index + node_offset

                all_edge_index.append(edge_index_global)
                all_edge_attr.append(edge_attr)
                all_edge_labels.append(edge_labels)
                graph_boundaries.append(graph_boundaries[-1] + n_edges)
                node_offset += len(ev["hits"])

            if all_edge_index:
                save_dict["edge_index"] = np.concatenate(all_edge_index, axis=1)
                save_dict["edge_attr"] = np.concatenate(all_edge_attr, axis=0)
                save_dict["edge_labels"] = np.concatenate(all_edge_labels)
            else:
                save_dict["edge_index"] = np.zeros((2, 0), dtype=np.int64)
                save_dict["edge_attr"] = np.zeros((0, 5), dtype=np.float32)
                save_dict["edge_labels"] = np.zeros(0, dtype=np.float32)
            save_dict["graph_boundaries"] = np.array(graph_boundaries, dtype=np.int64)

        np.savez_compressed(output_path, **save_dict)

    elif fmt == "pt":
        import torch
        torch.save(events, output_path)


def main():
    args = parse_args()
    t_start = time.time()

    print("=" * 70)
    print("G4SBS Training Data Generator for SBS-offline AIML")
    print("=" * 70)

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load events ---
    events = []

    if args.synthetic:
        print(f"\nGenerating {args.n_events} synthetic events...")
        events = generate_synthetic_gem_data(
            n_events=args.n_events,
            n_layers=args.n_layers,
            n_signal_tracks=args.n_signal_tracks,
            n_noise_hits_per_layer=args.n_noise_hits,
            rng_seed=args.seed,
        )
        print(f"  Generated {len(events)} events")
    elif args.input:
        files = expand_input_files(args.input)
        if not files:
            print("Error: No input files found.")
            sys.exit(1)

        print(f"\nLoading G4SBS data from {len(files)} file(s)...")
        for fpath in files:
            print(f"  Reading: {fpath}")
            try:
                file_events, mc_tracks = load_g4sbs_gem_data(
                    fpath,
                    tree_name=args.tree_name,
                    det_prefix=args.det_prefix,
                    max_events=args.max_events,
                )
                events.extend(file_events)
                print(f"    -> {len(file_events)} events loaded")
            except Exception as e:
                print(f"    -> Error: {e}")
                continue

        if not events:
            print("Error: No events loaded from any file.")
            sys.exit(1)
    else:
        print("Error: Specify --input or --synthetic")
        sys.exit(1)

    # --- Re-label if needed ---
    if not args.synthetic and (args.min_hits_per_track != 3 or args.primary_only):
        print(f"\nRe-labeling with min_hits={args.min_hits_per_track}, "
              f"primary_only={args.primary_only}...")
        events = relabel_events(
            events,
            min_hits_per_track=args.min_hits_per_track,
            primary_only=args.primary_only,
        )

    # --- Compute statistics ---
    print("\nComputing dataset statistics...")
    stats, feat_mean, feat_std = compute_statistics(events)

    print(f"  Events:            {stats['n_events']}")
    print(f"  Total hits:        {stats['n_total_hits']}")
    print(f"  Signal hits:       {stats['n_signal_hits']} "
          f"({100*stats['signal_fraction']:.1f}%)")
    print(f"  Noise hits:        {stats['n_noise_hits']}")
    print(f"  Hits/event:        {stats['hits_per_event_mean']:.1f} "
          f"+/- {stats['hits_per_event_std']:.1f}")
    print(f"  Tracks/event:      {stats['tracks_per_event_mean']:.1f} "
          f"+/- {stats['tracks_per_event_std']:.1f}")
    print(f"  Layers:            {stats['n_layers']} {stats['layer_ids']}")

    # --- Load waveforms if requested ---
    waveform_events = None
    if args.include_waveforms and args.input:
        print("\nLoading digitized strip waveforms...")
        waveform_events = []
        files = expand_input_files(args.input)
        for fpath in files:
            try:
                wf_events = load_g4sbs_digitized_strips(
                    fpath,
                    tree_name=args.tree_name,
                    det_prefix=args.det_prefix,
                    max_events=args.max_events,
                )
                waveform_events.extend(wf_events)
                print(f"  {fpath}: {len(wf_events)} events, "
                      f"{sum(len(e['waveforms']) for e in wf_events)} strips")
            except Exception as e:
                print(f"  {fpath}: Error - {e}")
        if waveform_events:
            total_strips = sum(len(e["waveforms"]) for e in waveform_events)
            print(f"  Total waveforms: {total_strips}")

    # --- Split into train/val/test ---
    print(f"\nSplitting: train={args.train_fraction}, val={args.val_fraction}, "
          f"test={1 - args.train_fraction - args.val_fraction:.2f}")

    train_events, val_events, test_events = split_events(
        events, args.train_fraction, args.val_fraction, seed=args.seed
    )

    print(f"  Train: {len(train_events)} events")
    print(f"  Val:   {len(val_events)} events")
    print(f"  Test:  {len(test_events)} events")

    # --- Save ---
    ext = ".npz" if args.format == "npz" else ".pt"

    print(f"\nSaving to {args.output_dir}/...")
    for split_name, split_events_data in [
        ("train", train_events), ("val", val_events), ("test", test_events)
    ]:
        output_path = os.path.join(args.output_dir, f"{split_name}{ext}")
        save_split(
            split_events_data, output_path, fmt=args.format,
            build_graphs=args.build_graphs,
            r_connect=args.r_connect, k_nearest=args.k_nearest,
        )
        file_size = os.path.getsize(output_path) / 1e6
        print(f"  {split_name}: {output_path} ({file_size:.1f} MB)")

    # Save normalization parameters
    norm_path = os.path.join(args.output_dir, "normalization.npz")
    feat_std_safe = feat_std.copy()
    feat_std_safe[feat_std_safe < 1e-8] = 1.0
    np.savez(norm_path, mean=feat_mean, std=feat_std_safe)
    print(f"  normalization: {norm_path}")

    # Save statistics and config
    config = {
        "args": vars(args),
        "statistics": stats,
        "feature_names": [
            "x_mm", "y_mm", "z_mm", "edep_keV", "t_ns", "t_corr_ns",
            "slope_x", "slope_y", "momentum_GeV", "beta",
            "edep_deconv_keV", "slope_x_deconv", "t_deconv_ns", "t_rms_ns",
        ],
    }
    config_path = os.path.join(args.output_dir, "config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2, default=str)
    print(f"  config: {config_path}")

    # Save waveforms separately if loaded
    if waveform_events:
        wf_train, wf_val, wf_test = split_events(
            waveform_events, args.train_fraction, args.val_fraction, seed=args.seed
        )
        for split_name, wf_split in [
            ("train_waveforms", wf_train),
            ("val_waveforms", wf_val),
            ("test_waveforms", wf_test),
        ]:
            wf_path = os.path.join(args.output_dir, f"{split_name}{ext}")
            all_wf = np.concatenate([e["waveforms"] for e in wf_split if len(e["waveforms"]) > 0])
            np.savez_compressed(wf_path, waveforms=all_wf)
            print(f"  {split_name}: {wf_path} ({os.path.getsize(wf_path)/1e6:.1f} MB)")

    dt = time.time() - t_start
    print(f"\nDone in {dt:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
