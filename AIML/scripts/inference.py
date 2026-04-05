#!/usr/bin/env python3
"""
Inference pipeline for applying trained AIML models to SBS-offline data.

Runs the full pipeline:
1. Load GEM hits from ROOT file
2. Apply autoencoder denoising + hit classification
3. Build event graphs from cleaned hits
4. Apply GNN track finding
5. Export results back to ROOT format

Usage:
    python inference.py \
        --input /path/to/data.root \
        --autoencoder-checkpoint checkpoints/autoencoder/best_model.pt \
        --gnn-checkpoint checkpoints/gnn_tracker/best_model.pt \
        --output results.root

    # Run on synthetic data for testing:
    python inference.py --synthetic --output results_synthetic.npz
"""

import argparse
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.autoencoder import GEMHitAutoencoder
from models.gnn_tracker import TrackGNN, extract_tracks_from_edges
from data.root_io import (
    generate_synthetic_gem_data, hits_to_graph, HIT_FEATURE_NAMES,
)


def load_autoencoder(checkpoint_path, device):
    """Load trained autoencoder from checkpoint."""
    ckpt = torch.load(checkpoint_path, map_location=device)
    args = ckpt.get("args", {})

    model = GEMHitAutoencoder(
        n_features=args.get("n_features", 14),
        latent_dim=args.get("latent_dim", 8),
        hidden_dims=args.get("hidden_dims", [64, 32]),
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()

    # Load normalization
    norm_dir = os.path.dirname(checkpoint_path)
    norm_path = os.path.join(norm_dir, "normalization.npz")
    if os.path.exists(norm_path):
        norm = np.load(norm_path)
        return model, norm["mean"], norm["std"]

    return model, None, None


def load_gnn(checkpoint_path, device):
    """Load trained GNN tracker from checkpoint."""
    ckpt = torch.load(checkpoint_path, map_location=device)
    args = ckpt.get("args", {})

    model = TrackGNN(
        node_dim=args.get("node_dim", 14),
        edge_dim=args.get("edge_dim", 5),
        hidden_dim=args.get("hidden_dim", 64),
        n_iterations=args.get("n_iterations", 4),
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()

    norm_dir = os.path.dirname(checkpoint_path)
    norm_path = os.path.join(norm_dir, "normalization.npz")
    if os.path.exists(norm_path):
        norm = np.load(norm_path)
        return model, norm["mean"], norm["std"]

    return model, None, None


@torch.no_grad()
def denoise_hits(model, hits, mean, std, device, signal_threshold=0.5):
    """Apply autoencoder denoising and classification to hit features.

    Parameters
    ----------
    model : GEMHitAutoencoder
        Trained autoencoder model.
    hits : np.ndarray
        Shape (n_hits, n_features).
    mean, std : np.ndarray
        Normalization parameters.
    device : torch.device
    signal_threshold : float
        Classification threshold for signal/noise.

    Returns
    -------
    denoised_hits : np.ndarray
        Denoised feature vectors.
    signal_scores : np.ndarray
        Signal probability for each hit.
    signal_mask : np.ndarray
        Boolean mask of classified signal hits.
    """
    # Normalize
    if mean is not None:
        hits_norm = (hits - mean) / std
    else:
        hits_norm = hits

    x = torch.tensor(hits_norm, dtype=torch.float32).to(device)
    x_recon, z, scores = model(x)

    # Denormalize reconstructed hits
    denoised = x_recon.cpu().numpy()
    if mean is not None:
        denoised = denoised * std + mean

    signal_probs = torch.sigmoid(scores).cpu().numpy()
    signal_mask = signal_probs > signal_threshold

    return denoised, signal_probs, signal_mask


@torch.no_grad()
def find_tracks(model, hits, layer_ids, mean, std, device,
                r_connect=50.0, k_nearest=8, edge_threshold=0.5):
    """Apply GNN track finder to an event.

    Parameters
    ----------
    model : TrackGNN
        Trained GNN model.
    hits : np.ndarray
        Shape (n_hits, n_features).
    layer_ids : np.ndarray
        Layer assignments for each hit.
    mean, std : np.ndarray
        Normalization parameters.
    device : torch.device

    Returns
    -------
    track_candidates : list[list[int]]
        Track candidates as lists of hit indices.
    edge_scores : np.ndarray
        Edge classification scores.
    node_scores : np.ndarray
        Node signal/noise scores.
    """
    # Build graph
    edge_index, edge_attr = hits_to_graph(
        hits, layer_ids, r_connect=r_connect, k_nearest=k_nearest
    )

    if edge_index.shape[1] == 0:
        return [], np.array([]), np.array([])

    # Normalize node features
    if mean is not None:
        hits_norm = (hits - mean) / std
    else:
        hits_norm = hits

    node_feat = torch.tensor(hits_norm, dtype=torch.float32).to(device)
    edge_idx = torch.tensor(edge_index, dtype=torch.long).to(device)
    edge_at = torch.tensor(edge_attr, dtype=torch.float32).to(device)

    edge_scores, node_scores, _ = model(node_feat, edge_idx, edge_at)
    edge_probs = torch.sigmoid(edge_scores)
    node_probs = torch.sigmoid(node_scores)

    # Extract tracks
    track_candidates = extract_tracks_from_edges(
        edge_idx, edge_probs, threshold=edge_threshold
    )

    return track_candidates, edge_probs.cpu().numpy(), node_probs.cpu().numpy()


def run_full_pipeline(events, autoencoder, ae_mean, ae_std,
                      gnn, gnn_mean, gnn_std, device, config):
    """Run the complete denoising + tracking pipeline on a list of events.

    Parameters
    ----------
    events : list[dict]
        Event data with "hits", "layer_ids", and optionally "labels", "track_ids".
    autoencoder : GEMHitAutoencoder or None
    gnn : TrackGNN or None
    device : torch.device
    config : dict
        Pipeline configuration.

    Returns
    -------
    results : list[dict]
        Per-event results with denoised hits, signal scores, and tracks.
    """
    results = []

    for i, event in enumerate(events):
        hits = event["hits"]
        layer_ids = event["layer_ids"]
        result = {"event_idx": i, "n_input_hits": len(hits)}

        # Step 1: Autoencoder denoising
        if autoencoder is not None:
            denoised, signal_scores, signal_mask = denoise_hits(
                autoencoder, hits, ae_mean, ae_std, device,
                signal_threshold=config.get("signal_threshold", 0.5),
            )
            result["denoised_hits"] = denoised
            result["signal_scores"] = signal_scores
            result["signal_mask"] = signal_mask
            result["n_signal_hits"] = signal_mask.sum()

            # Use denoised signal hits for tracking
            clean_hits = denoised[signal_mask]
            clean_layers = layer_ids[signal_mask]
        else:
            clean_hits = hits
            clean_layers = layer_ids

        # Step 2: GNN track finding
        if gnn is not None and len(clean_hits) >= 3:
            tracks, edge_scores, node_scores = find_tracks(
                gnn, clean_hits, clean_layers, gnn_mean, gnn_std, device,
                r_connect=config.get("r_connect", 50.0),
                k_nearest=config.get("k_nearest", 8),
                edge_threshold=config.get("edge_threshold", 0.5),
            )
            result["tracks"] = tracks
            result["n_tracks"] = len(tracks)
            result["edge_scores"] = edge_scores
            result["node_scores"] = node_scores
        else:
            result["tracks"] = []
            result["n_tracks"] = 0

        results.append(result)

    return results


def parse_args():
    parser = argparse.ArgumentParser(description="AIML inference pipeline")
    parser.add_argument("--input", type=str, default=None,
                        help="Input ROOT file")
    parser.add_argument("--synthetic", action="store_true",
                        help="Run on synthetic data")
    parser.add_argument("--n-events", type=int, default=100)
    parser.add_argument("--autoencoder-checkpoint", type=str, default=None,
                        help="Autoencoder checkpoint path")
    parser.add_argument("--gnn-checkpoint", type=str, default=None,
                        help="GNN checkpoint path")
    parser.add_argument("--signal-threshold", type=float, default=0.5)
    parser.add_argument("--edge-threshold", type=float, default=0.5)
    parser.add_argument("--r-connect", type=float, default=50.0)
    parser.add_argument("--k-nearest", type=int, default=8)
    parser.add_argument("--output", type=str, default="aiml_results.npz",
                        help="Output file path (.npz or .root)")
    parser.add_argument("--device", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load models
    autoencoder, ae_mean, ae_std = None, None, None
    gnn, gnn_mean, gnn_std = None, None, None

    if args.autoencoder_checkpoint:
        print(f"Loading autoencoder from {args.autoencoder_checkpoint}")
        autoencoder, ae_mean, ae_std = load_autoencoder(
            args.autoencoder_checkpoint, device
        )

    if args.gnn_checkpoint:
        print(f"Loading GNN from {args.gnn_checkpoint}")
        gnn, gnn_mean, gnn_std = load_gnn(args.gnn_checkpoint, device)

    if autoencoder is None and gnn is None:
        print("Warning: No model checkpoints provided. Running in passthrough mode.")

    # Load data
    if args.synthetic:
        print(f"Generating {args.n_events} synthetic events...")
        events = generate_synthetic_gem_data(n_events=args.n_events)
    else:
        print("ROOT file input for inference pipeline requires uproot.")
        sys.exit(1)

    # Run pipeline
    config = {
        "signal_threshold": args.signal_threshold,
        "edge_threshold": args.edge_threshold,
        "r_connect": args.r_connect,
        "k_nearest": args.k_nearest,
    }

    print(f"Running inference on {len(events)} events...")
    t0 = time.time()
    results = run_full_pipeline(
        events, autoencoder, ae_mean, ae_std, gnn, gnn_mean, gnn_std, device, config
    )
    dt = time.time() - t0
    print(f"Inference complete in {dt:.1f}s ({len(events)/dt:.0f} events/s)")

    # Summary statistics
    n_tracks_total = sum(r["n_tracks"] for r in results)
    print(f"Total tracks found: {n_tracks_total}")
    print(f"Average tracks/event: {n_tracks_total/len(results):.1f}")

    if "signal_mask" in results[0]:
        n_signal = sum(r["n_signal_hits"] for r in results)
        n_total = sum(r["n_input_hits"] for r in results)
        print(f"Signal hits: {n_signal}/{n_total} ({100*n_signal/n_total:.1f}%)")

    # Save results
    if args.output.endswith(".npz"):
        save_data = {
            "n_events": len(results),
            "n_tracks": np.array([r["n_tracks"] for r in results]),
            "n_input_hits": np.array([r["n_input_hits"] for r in results]),
        }
        np.savez(args.output, **save_data)
        print(f"Results saved to {args.output}")
    else:
        print(f"Results summary printed above. ROOT output not yet implemented.")


if __name__ == "__main__":
    main()
