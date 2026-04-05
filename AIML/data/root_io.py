"""
ROOT I/O utilities for extracting GEM hit and track data from SBS-offline
ROOT files and converting to PyTorch-compatible formats.

Supports both real data and simulation (G4SBS) ROOT tree formats.
"""

import numpy as np

try:
    import uproot
except ImportError:
    uproot = None

try:
    import ROOT
    _HAS_PYROOT = True
except ImportError:
    _HAS_PYROOT = False


# ---------------------------------------------------------------------------
# Feature definitions matching SBS-offline GEM data structures
# ---------------------------------------------------------------------------

# Per-hit features extracted from sbsgemhit_t
HIT_FEATURE_NAMES = [
    "xghit",       # global X coordinate
    "yghit",       # global Y coordinate
    "zghit",       # global Z coordinate
    "Ehit",        # ADC energy (average of U/V cluster sums)
    "thit",        # ADC-weighted mean time
    "thitcorr",    # corrected hit time
    "ADCasym",     # (ADCU - ADCV) / (ADCU + ADCV)
    "tdiff",       # time difference tu - tv
    "corrcoeff_clust",   # cluster-level XY correlation
    "corrcoeff_strip",   # strip-level XY correlation
    "EhitDeconv",        # deconvoluted energy
    "ADCasymDeconv",     # deconvoluted ADC asymmetry
    "thitDeconv",        # deconvoluted hit time
    "tdiffDeconv",       # deconvoluted time difference
]

# Per-cluster (1D) features from sbsgemcluster_t
CLUSTER_FEATURE_NAMES = [
    "nstrips",
    "hitpos_mean",
    "hitpos_sigma",
    "clusterADCsum",
    "clusterADCsumDeconv",
    "t_mean",
    "t_mean_deconv",
]


def _check_uproot():
    if uproot is None:
        raise ImportError(
            "uproot is required for ROOT I/O. Install with: pip install uproot"
        )


def load_gem_hits_from_root(
    filepath,
    tree_name="T",
    branch_prefix="sbs.gems",
    max_events=None,
    features=None,
):
    """Load GEM 2D hit data from a SBS-offline ROOT file.

    Parameters
    ----------
    filepath : str
        Path to the ROOT file.
    tree_name : str
        Name of the TTree (default: "T").
    branch_prefix : str
        Branch prefix for GEM hit arrays (default: "sbs.gems").
    max_events : int or None
        Maximum number of events to read.
    features : list[str] or None
        Subset of HIT_FEATURE_NAMES to load. None loads all.

    Returns
    -------
    list[np.ndarray]
        List of per-event arrays, each with shape (n_hits, n_features).
    list[np.ndarray] or None
        List of per-event truth labels (1=signal, 0=noise) if MC truth
        branches are available, otherwise None.
    """
    _check_uproot()

    if features is None:
        features = HIT_FEATURE_NAMES

    branch_names = [f"{branch_prefix}.hit.{f}" for f in features]

    # Also try to load MC truth if available
    mc_branch = f"{branch_prefix}.hit.ontrack"

    with uproot.open(filepath) as f:
        tree = f[tree_name]

        available = set(tree.keys())
        load_branches = [b for b in branch_names if b in available]

        if not load_branches:
            raise ValueError(
                f"No matching branches found. Available: {sorted(available)[:20]}..."
            )

        has_mc = mc_branch in available
        if has_mc:
            load_branches.append(mc_branch)

        entry_stop = max_events if max_events else None
        arrays = tree.arrays(load_branches, entry_stop=entry_stop, library="np")

        events_hits = []
        events_labels = []
        n_events = len(next(iter(arrays.values())))

        for i in range(n_events):
            hit_data = []
            for b in branch_names:
                if b in arrays:
                    hit_data.append(np.asarray(arrays[b][i], dtype=np.float32))
            if hit_data:
                event_array = np.stack(hit_data, axis=-1)
                events_hits.append(event_array)
                if has_mc:
                    events_labels.append(
                        np.asarray(arrays[mc_branch][i], dtype=np.float32)
                    )

    labels = events_labels if has_mc else None
    return events_hits, labels


def load_tracks_from_root(
    filepath,
    tree_name="T",
    branch_prefix="sbs.gems",
    max_events=None,
):
    """Load reconstructed track parameters from a SBS-offline ROOT file.

    Returns per-event arrays of track parameters (x, y, xp, yp, chi2, ndf).
    """
    _check_uproot()

    track_branches = [
        f"{branch_prefix}.track.x",
        f"{branch_prefix}.track.y",
        f"{branch_prefix}.track.xp",
        f"{branch_prefix}.track.yp",
        f"{branch_prefix}.track.chi2ndf",
        f"{branch_prefix}.track.nhits",
    ]

    with uproot.open(filepath) as f:
        tree = f[tree_name]
        available = set(tree.keys())
        load_branches = [b for b in track_branches if b in available]

        if not load_branches:
            return []

        entry_stop = max_events if max_events else None
        arrays = tree.arrays(load_branches, entry_stop=entry_stop, library="np")

        n_events = len(next(iter(arrays.values())))
        events_tracks = []

        for i in range(n_events):
            track_data = []
            for b in load_branches:
                track_data.append(np.asarray(arrays[b][i], dtype=np.float32))
            if track_data:
                events_tracks.append(np.stack(track_data, axis=-1))

    return events_tracks


def hits_to_graph(hits, layer_ids, r_connect=50.0, k_nearest=8):
    """Convert per-event hit arrays to a graph representation for GNN tracking.

    Builds a graph where:
    - Nodes are detector hits (with hit features as node features)
    - Edges connect hits in adjacent detector layers that are
      geometrically compatible (within r_connect or k-nearest neighbors)

    Parameters
    ----------
    hits : np.ndarray
        Shape (n_hits, n_features). First 3 features must be (x, y, z).
    layer_ids : np.ndarray
        Integer layer ID for each hit, shape (n_hits,).
    r_connect : float
        Maximum transverse distance (mm) to form an edge between hits
        in adjacent layers.
    k_nearest : int
        Maximum number of nearest-neighbor edges per hit.

    Returns
    -------
    edge_index : np.ndarray
        Shape (2, n_edges) - source and target node indices.
    edge_attr : np.ndarray
        Shape (n_edges, n_edge_features) - edge features (dx, dy, dz, dr, dphi).
    """
    n_hits = len(hits)
    if n_hits == 0:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 5), dtype=np.float32)

    positions = hits[:, :3]  # x, y, z
    unique_layers = np.unique(layer_ids)
    sorted_layers = np.sort(unique_layers)

    src_list = []
    dst_list = []
    edge_feat_list = []

    # Connect hits between adjacent layers
    for i in range(len(sorted_layers) - 1):
        layer_a = sorted_layers[i]
        layer_b = sorted_layers[i + 1]

        mask_a = layer_ids == layer_a
        mask_b = layer_ids == layer_b

        idx_a = np.where(mask_a)[0]
        idx_b = np.where(mask_b)[0]

        if len(idx_a) == 0 or len(idx_b) == 0:
            continue

        pos_a = positions[idx_a]
        pos_b = positions[idx_b]

        # Compute pairwise transverse distances (XY plane)
        dx = pos_b[:, 0][np.newaxis, :] - pos_a[:, 0][:, np.newaxis]
        dy = pos_b[:, 1][np.newaxis, :] - pos_a[:, 1][:, np.newaxis]
        dz = pos_b[:, 2][np.newaxis, :] - pos_a[:, 2][:, np.newaxis]
        dr = np.sqrt(dx ** 2 + dy ** 2)

        for ia in range(len(idx_a)):
            distances = dr[ia]
            # Apply radius cut
            within_r = distances < r_connect
            candidate_indices = np.where(within_r)[0]

            if len(candidate_indices) == 0:
                continue

            # Keep at most k_nearest
            if len(candidate_indices) > k_nearest:
                sorted_idx = np.argsort(distances[candidate_indices])
                candidate_indices = candidate_indices[sorted_idx[:k_nearest]]

            for ib in candidate_indices:
                global_a = idx_a[ia]
                global_b = idx_b[ib]

                ddx = dx[ia, ib]
                ddy = dy[ia, ib]
                ddz = dz[ia, ib]
                ddr = dr[ia, ib]
                dphi = np.arctan2(ddy, ddx)

                # Bidirectional edges
                src_list.extend([global_a, global_b])
                dst_list.extend([global_b, global_a])
                feat = [ddx, ddy, ddz, ddr, dphi]
                edge_feat_list.append(feat)
                edge_feat_list.append([-ddx, -ddy, -ddz, ddr, dphi + np.pi])

    if len(src_list) == 0:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 5), dtype=np.float32)

    edge_index = np.array([src_list, dst_list], dtype=np.int64)
    edge_attr = np.array(edge_feat_list, dtype=np.float32)

    return edge_index, edge_attr


# ---------------------------------------------------------------------------
# G4SBS Monte Carlo truth tree formats
# ---------------------------------------------------------------------------

# Branch naming conventions for G4SBS truth-level GEM hit trees.
# The tree name is "T" and branches follow the pattern:
#   {detector_prefix}_hit_{field}
# Common detector prefixes for GEM:
#   Earm.BBGEM   - BigBite GEM (GMn, GEn experiments)
#
# These provide MC truth needed for supervised training labels.

G4SBS_GEM_HIT_FIELDS = [
    "plane",   # GEM plane/layer index
    "x",       # local x position (m)
    "y",       # local y position (m)
    "z",       # local z position (m)
    "xg",      # global x position (m)
    "yg",      # global y position (m)
    "zg",      # global z position (m)
    "t",       # average hit time (ns)
    "trms",    # hit time RMS (ns)
    "tmin",    # minimum hit time (ns)
    "tmax",    # maximum hit time (ns)
    "tx",      # local track x at plane (m)
    "ty",      # local track y at plane (m)
    "txp",     # track dx/dz slope
    "typ",     # track dy/dz slope
    "trid",    # Geant4 track ID
    "mid",     # Geant4 mother track ID
    "pid",     # PDG particle ID
    "p",       # particle momentum (GeV)
    "edep",    # energy deposit (GeV)
    "beta",    # particle velocity v/c
    "vx",      # vertex x (m)
    "vy",      # vertex y (m)
    "vz",      # vertex z (m)
    "xin",     # entry point x (m)
    "yin",     # entry point y (m)
    "zin",     # entry point z (m)
    "xout",    # exit point x (m)
    "yout",    # exit point y (m)
    "zout",    # exit point z (m)
]

G4SBS_GEM_TRACK_FIELDS = [
    "TID",       # Geant4 track ID
    "PID",       # PDG particle ID
    "MID",       # mother track ID
    "NumHits",   # number of hits on track
    "NumPlanes", # number of planes hit
    "NDF",       # degrees of freedom
    "Chi2fit",   # chi2 of fit
    "Chi2true",  # chi2 vs true trajectory
    "X",         # track x at reference plane (m)
    "Y",         # track y at reference plane (m)
    "Xp",        # track dx/dz
    "Yp",        # track dy/dz
    "T",         # track time (ns)
    "P",         # track momentum (GeV)
    "Sx",        # spin x component
    "Sy",        # spin y component
    "Sz",        # spin z component
    "Xfit",      # fitted track x (m)
    "Yfit",      # fitted track y (m)
    "Xpfit",     # fitted dx/dz
    "Ypfit",     # fitted dy/dz
]

# Digitized GEM strip branches (6 ADC time samples per strip)
G4SBS_GEM_DIG_FIELDS = ["strip", "adc_0", "adc_1", "adc_2", "adc_3", "adc_4", "adc_5"]


def load_g4sbs_gem_data(
    filepath,
    tree_name="T",
    det_prefix="Earm_BBGEM",
    max_events=None,
    convert_to_mm=True,
):
    """Load G4SBS truth-level GEM hit and track data for training.

    Reads the Geant4 simulation truth tree produced by G4SBS and converts
    it into the per-event dict format expected by the AIML training pipeline.

    Parameters
    ----------
    filepath : str
        Path to G4SBS output ROOT file (or digitized file).
    tree_name : str
        TTree name (default "T").
    det_prefix : str
        Detector branch prefix. Use "Earm_BBGEM" for BigBite GEM (GMn/GEn).
        Branch names follow the pattern ``{det_prefix}_hit_{field}``.
        Note: dots in branch names may appear as underscores in uproot.
    max_events : int or None
        Maximum number of events to load.
    convert_to_mm : bool
        Convert G4SBS positions from meters to mm (the SBS-offline convention).

    Returns
    -------
    events : list[dict]
        Per-event dictionaries compatible with GEMTrackGraphDataset:
        - "hits": (n_hits, n_features) float32 array
        - "layer_ids": (n_hits,) int64 array of plane indices
        - "labels": (n_hits,) float32 (1=primary signal, 0=secondary/noise)
        - "track_ids": (n_hits,) int64 Geant4 track IDs
        - "true_tracks": list of dicts with track parameters
        - "pids": (n_hits,) int64 PDG particle IDs
        - "momenta": (n_hits,) float32 particle momenta (GeV)
        - "edep": (n_hits,) float32 energy deposits (GeV)
    mc_tracks : list[list[dict]]
        Per-event list of MC truth track parameter dicts.
    """
    _check_uproot()

    # Build branch names - try both dot and underscore separators
    hit_prefix = f"{det_prefix}_hit"
    track_prefix = f"{det_prefix}_Track"
    nhits_branch = f"{det_prefix}_hit_nhits"
    ntracks_branch = f"{det_prefix}_Track_ntracks"

    hit_branches = [f"{hit_prefix}_{f}" for f in G4SBS_GEM_HIT_FIELDS]
    track_branches = [f"{track_prefix}_{f}" for f in G4SBS_GEM_TRACK_FIELDS]

    with uproot.open(filepath) as f:
        tree = f[tree_name]
        available = set(tree.keys())

        # Auto-detect branch naming: dots vs underscores
        if nhits_branch not in available:
            alt_prefix = det_prefix.replace("_", ".")
            hit_prefix_alt = f"{alt_prefix}.hit"
            track_prefix_alt = f"{alt_prefix}.Track"
            nhits_branch_alt = f"{alt_prefix}.hit.nhits"

            if nhits_branch_alt in available:
                hit_branches = [f"{hit_prefix_alt}.{f}" for f in G4SBS_GEM_HIT_FIELDS]
                track_branches = [f"{track_prefix_alt}.{f}" for f in G4SBS_GEM_TRACK_FIELDS]
                nhits_branch = nhits_branch_alt
                ntracks_branch = f"{alt_prefix}.Track.ntracks"

        # Filter to available branches
        hit_branches_avail = [b for b in hit_branches if b in available]
        track_branches_avail = [b for b in track_branches if b in available]

        if not hit_branches_avail:
            raise ValueError(
                f"No G4SBS GEM hit branches found with prefix '{det_prefix}'. "
                f"Available branches (first 30): {sorted(available)[:30]}"
            )

        all_branches = hit_branches_avail[:]
        if nhits_branch in available:
            all_branches.append(nhits_branch)
        all_branches.extend(track_branches_avail)
        if ntracks_branch in available:
            all_branches.append(ntracks_branch)

        entry_stop = max_events if max_events else None
        arrays = tree.arrays(all_branches, entry_stop=entry_stop, library="np")

        n_events = len(next(iter(arrays.values())))
        scale = 1000.0 if convert_to_mm else 1.0  # m -> mm

        # Build a field lookup for hit data
        def _get_hit(field, event_idx):
            for sep in ["_", "."]:
                key = f"{hit_prefix}{sep}{field}" if sep == "_" else f"{hit_prefix.replace('_', '.')}.{field}"
                # Try both naming conventions
                for candidate in [f"{hit_prefix}_{field}",
                                  f"{hit_prefix.replace('_hit', '.hit')}.{field}"]:
                    if candidate in arrays:
                        return np.asarray(arrays[candidate][event_idx])
            return None

        def _get_track(field, event_idx):
            for candidate in [f"{track_prefix}_{field}",
                              f"{track_prefix.replace('_Track', '.Track')}.{field}"]:
                if candidate in arrays:
                    return np.asarray(arrays[candidate][event_idx])
            return None

        events = []
        mc_tracks_all = []

        for ev in range(n_events):
            # --- Extract hit data ---
            plane = _get_hit("plane", ev)
            if plane is None or len(plane) == 0:
                continue

            n_hits = len(plane)

            # Positions (prefer global, fall back to local)
            xg = _get_hit("xg", ev)
            yg = _get_hit("yg", ev)
            zg = _get_hit("zg", ev)
            if xg is None:
                xg = _get_hit("x", ev)
                yg = _get_hit("y", ev)
                zg = _get_hit("z", ev)

            # Timing
            t_hit = _get_hit("t", ev)
            trms = _get_hit("trms", ev)

            # Track slopes
            txp = _get_hit("txp", ev)
            typ = _get_hit("typ", ev)

            # Track IDs and particle info
            trid = _get_hit("trid", ev)
            pid = _get_hit("pid", ev)
            momentum = _get_hit("p", ev)
            edep = _get_hit("edep", ev)
            beta = _get_hit("beta", ev)

            # Build feature matrix matching our 14-feature convention:
            # [x, y, z, E, t, t_corr, asym, tdiff, corr_clust, corr_strip,
            #  E_deconv, asym_deconv, t_deconv, tdiff_deconv]
            # For G4SBS truth data we fill what we have and use placeholders
            # for features that only exist after digitization + reconstruction.
            x_mm = xg * scale if xg is not None else np.zeros(n_hits)
            y_mm = yg * scale if yg is not None else np.zeros(n_hits)
            z_mm = zg * scale if zg is not None else np.zeros(n_hits)
            e_gev = edep if edep is not None else np.zeros(n_hits)
            t_ns = t_hit if t_hit is not None else np.zeros(n_hits)
            t_rms = trms if trms is not None else np.zeros(n_hits)
            slope_x = txp if txp is not None else np.zeros(n_hits)
            slope_y = typ if typ is not None else np.zeros(n_hits)
            mom = momentum if momentum is not None else np.zeros(n_hits)
            b = beta if beta is not None else np.zeros(n_hits)

            hits = np.column_stack([
                x_mm, y_mm, z_mm,          # positions (mm)
                e_gev * 1e6,                # energy deposit (keV, more natural scale)
                t_ns,                        # hit time (ns)
                t_ns,                        # corrected time (same as t for truth)
                slope_x,                     # dx/dz slope (proxy for ADC asymmetry)
                slope_y,                     # dy/dz slope (proxy for time diff)
                mom,                         # momentum (proxy for cluster correlation)
                b,                           # beta (proxy for strip correlation)
                e_gev * 1e6 * 0.9,          # placeholder for deconv energy
                slope_x,                     # placeholder for deconv asymmetry
                t_ns,                        # placeholder for deconv time
                t_rms,                       # time RMS as proxy for deconv tdiff
            ]).astype(np.float32)

            layer_ids = plane.astype(np.int64)

            # --- Truth labels ---
            # Primary signal: track ID 1 is typically the primary scattered particle
            # in G4SBS. Hits with trid > 0 and pid == +/-11 (electron) or
            # pid matching the primary are signal.
            trid_arr = trid.astype(np.int64) if trid is not None else -np.ones(n_hits, dtype=np.int64)
            pid_arr = pid.astype(np.int64) if pid is not None else np.zeros(n_hits, dtype=np.int64)

            # Signal = any hit from a particle that traverses >= 3 planes
            # (mimics the SBS-offline tracking requirement)
            unique_trids, trid_counts = np.unique(trid_arr[trid_arr > 0], return_counts=True)
            signal_trids = set(unique_trids[trid_counts >= 3])
            labels = np.array([1.0 if int(tid) in signal_trids else 0.0
                               for tid in trid_arr], dtype=np.float32)

            # --- Extract MC truth tracks ---
            mc_tracks = []
            tid_track = _get_track("TID", ev)
            if tid_track is not None and len(tid_track) > 0:
                for it in range(len(tid_track)):
                    track_dict = {}
                    for field in G4SBS_GEM_TRACK_FIELDS:
                        val = _get_track(field, ev)
                        if val is not None and it < len(val):
                            track_dict[field] = float(val[it])
                    # Convert position fields to mm
                    for pos_field in ["X", "Y", "Xfit", "Yfit"]:
                        if pos_field in track_dict:
                            track_dict[pos_field] *= scale
                    mc_tracks.append(track_dict)

            events.append({
                "hits": hits,
                "layer_ids": layer_ids,
                "labels": labels,
                "track_ids": trid_arr,
                "true_tracks": mc_tracks,
                "pids": pid_arr,
                "momenta": mom.astype(np.float32),
                "edep": (e_gev * 1e6).astype(np.float32),
            })
            mc_tracks_all.append(mc_tracks)

    return events, mc_tracks_all


def load_g4sbs_digitized_strips(
    filepath,
    tree_name="T",
    det_prefix="Earm_BBGEM",
    n_modules=5,
    max_events=None,
):
    """Load digitized GEM strip ADC samples from a G4SBS digitized file.

    Reads the 6-sample ADC waveforms per strip, organized by GEM module.
    Useful for training the GEMConvAutoencoder on waveform denoising.

    Parameters
    ----------
    filepath : str
        Path to digitized ROOT file.
    tree_name : str
        TTree name.
    det_prefix : str
        Detector prefix (e.g. "Earm_BBGEM").
    n_modules : int
        Number of GEM modules (planes * 2 for x/y readout).
        For BigBite GEM with 5 planes: modules 1x,1y,...,5x,5y.
    max_events : int or None
        Maximum events to read.

    Returns
    -------
    list[dict]
        Per-event dicts with:
        - "strips": list of (strip_id, module_id, adc_samples[6]) arrays per module
        - "waveforms": (n_strips, 6) float32 array of all ADC waveforms
        - "module_ids": (n_strips,) int array of module assignments
    """
    _check_uproot()

    with uproot.open(filepath) as f:
        tree = f[tree_name]
        available = set(tree.keys())

        # Detect branch naming for digitized GEM strips
        # Pattern: {prefix}_{layer}{axis}_dighit_{field}
        # e.g. Earm_BBGEM_1x_dighit_strip, Earm_BBGEM_1x_dighit_adc_0
        # or:  Earm.BBGEM.dighit.{field} (flat module-based)

        module_branches = {}
        for layer in range(1, 6):  # layers 1-5
            for axis in ["x", "y"]:
                mod_name = f"{layer}{axis}"
                mod_prefix = f"{det_prefix}_{mod_name}_dighit"
                # Also try dot notation
                mod_prefix_dot = f"{det_prefix.replace('_', '.')}_{mod_name}_dighit"

                for prefix_try in [mod_prefix, mod_prefix_dot]:
                    strip_key = f"{prefix_try}_strip"
                    if strip_key in available:
                        branches = [f"{prefix_try}_{f}" for f in G4SBS_GEM_DIG_FIELDS]
                        branches_avail = [b for b in branches if b in available]
                        if branches_avail:
                            module_branches[mod_name] = branches_avail
                        break

        # Also check flat dighit format (module + strip + adc_0..5)
        flat_prefix = f"{det_prefix}_dighit"
        flat_prefix_dot = f"{det_prefix.replace('_', '.')}.dighit"
        for prefix_try in [flat_prefix, flat_prefix_dot]:
            mod_key = f"{prefix_try}_module"
            if mod_key in available:
                branches = [f"{prefix_try}_{f}" for f in
                            ["module"] + G4SBS_GEM_DIG_FIELDS]
                branches_avail = [b for b in branches if b in available]
                if branches_avail:
                    module_branches["flat"] = branches_avail
                break

        if not module_branches:
            raise ValueError(
                f"No digitized GEM strip branches found with prefix '{det_prefix}'. "
                f"Available (first 30): {sorted(available)[:30]}"
            )

        # Load all needed branches
        all_branches = []
        for branches in module_branches.values():
            all_branches.extend(branches)
        all_branches = list(set(all_branches))

        entry_stop = max_events if max_events else None
        arrays = tree.arrays(all_branches, entry_stop=entry_stop, library="np")
        n_events_total = len(next(iter(arrays.values())))

        events = []
        for ev in range(n_events_total):
            all_waveforms = []
            all_module_ids = []

            if "flat" in module_branches:
                # Flat format: single set of branches with module index
                prefix_try = flat_prefix if f"{flat_prefix}_module" in arrays else flat_prefix_dot
                modules = np.asarray(arrays[f"{prefix_try}_module"][ev])
                strips = np.asarray(arrays[f"{prefix_try}_strip"][ev])
                adcs = []
                for s in range(6):
                    key = f"{prefix_try}_adc_{s}"
                    if key in arrays:
                        adcs.append(np.asarray(arrays[key][ev], dtype=np.float32))
                    else:
                        adcs.append(np.zeros(len(modules), dtype=np.float32))
                waveforms = np.stack(adcs, axis=-1)
                all_waveforms.append(waveforms)
                all_module_ids.append(modules.astype(np.int32))
            else:
                # Per-module format
                mod_idx = 0
                for mod_name, branches in sorted(module_branches.items()):
                    strip_key = [b for b in branches if b.endswith("_strip")]
                    if not strip_key:
                        continue
                    strips = np.asarray(arrays[strip_key[0]][ev])
                    n_strips_ev = len(strips)
                    if n_strips_ev == 0:
                        mod_idx += 1
                        continue

                    adcs = []
                    for s in range(6):
                        adc_key = [b for b in branches if b.endswith(f"_adc_{s}")]
                        if adc_key and adc_key[0] in arrays:
                            adcs.append(np.asarray(arrays[adc_key[0]][ev], dtype=np.float32))
                        else:
                            adcs.append(np.zeros(n_strips_ev, dtype=np.float32))
                    waveforms = np.stack(adcs, axis=-1)
                    all_waveforms.append(waveforms)
                    all_module_ids.append(np.full(n_strips_ev, mod_idx, dtype=np.int32))
                    mod_idx += 1

            if all_waveforms:
                events.append({
                    "waveforms": np.concatenate(all_waveforms, axis=0),
                    "module_ids": np.concatenate(all_module_ids, axis=0),
                })
            else:
                events.append({"waveforms": np.zeros((0, 6), dtype=np.float32),
                               "module_ids": np.zeros(0, dtype=np.int32)})

    return events


def generate_synthetic_gem_data(
    n_events=1000,
    n_layers=5,
    n_signal_tracks=2,
    n_noise_hits_per_layer=10,
    layer_spacing=100.0,
    hit_sigma=2.0,
    rng_seed=42,
):
    """Generate synthetic GEM-like hit data for development and testing.

    Creates straight-line tracks through multiple layers with Gaussian
    smearing and uniform random noise hits.

    Parameters
    ----------
    n_events : int
        Number of events to generate.
    n_layers : int
        Number of detector layers.
    n_signal_tracks : int
        Number of true tracks per event.
    n_noise_hits_per_layer : int
        Number of random noise hits per layer.
    layer_spacing : float
        Z spacing between layers (mm).
    hit_sigma : float
        Gaussian smearing of hit positions (mm).
    rng_seed : int
        Random seed for reproducibility.

    Returns
    -------
    list[dict]
        Per-event dictionaries with keys:
        - "hits": (n_hits, n_features) array
        - "layer_ids": (n_hits,) integer layer assignments
        - "labels": (n_hits,) 1=signal, 0=noise
        - "track_ids": (n_hits,) track index (-1 for noise)
        - "true_tracks": list of (x0, y0, dxdz, dydz) tuples
    """
    rng = np.random.RandomState(rng_seed)
    events = []

    z_layers = np.arange(n_layers) * layer_spacing
    detector_half_size = 300.0  # mm

    for _ in range(n_events):
        all_hits = []
        all_layers = []
        all_labels = []
        all_track_ids = []
        true_tracks = []

        # Generate signal tracks
        for t in range(n_signal_tracks):
            x0 = rng.uniform(-100, 100)
            y0 = rng.uniform(-100, 100)
            dxdz = rng.uniform(-0.1, 0.1)
            dydz = rng.uniform(-0.1, 0.1)
            true_tracks.append((x0, y0, dxdz, dydz))

            for layer, z in enumerate(z_layers):
                x = x0 + dxdz * z + rng.normal(0, hit_sigma)
                y = y0 + dydz * z + rng.normal(0, hit_sigma)

                # Simulated hit features: x, y, z, E, t, tcorr, asym, tdiff,
                # corr_clust, corr_strip, E_deconv, asym_deconv, t_deconv, tdiff_deconv
                energy = rng.exponential(500) + 200
                time = rng.normal(0, 5)
                hit = [
                    x, y, z,
                    energy, time, time + rng.normal(0, 1),
                    rng.normal(0, 0.1), rng.normal(0, 2),
                    rng.uniform(0.5, 1.0), rng.uniform(0.5, 1.0),
                    energy * 0.9, rng.normal(0, 0.1),
                    time + rng.normal(0, 0.5), rng.normal(0, 1),
                ]
                all_hits.append(hit)
                all_layers.append(layer)
                all_labels.append(1)
                all_track_ids.append(t)

        # Generate noise hits
        for layer, z in enumerate(z_layers):
            for _ in range(rng.poisson(n_noise_hits_per_layer)):
                x = rng.uniform(-detector_half_size, detector_half_size)
                y = rng.uniform(-detector_half_size, detector_half_size)
                energy = rng.exponential(100)
                time = rng.uniform(-50, 50)
                hit = [
                    x, y, z,
                    energy, time, time + rng.normal(0, 5),
                    rng.uniform(-1, 1), rng.uniform(-20, 20),
                    rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5),
                    energy * 0.8, rng.uniform(-1, 1),
                    time + rng.normal(0, 5), rng.uniform(-10, 10),
                ]
                all_hits.append(hit)
                all_layers.append(layer)
                all_labels.append(0)
                all_track_ids.append(-1)

        events.append({
            "hits": np.array(all_hits, dtype=np.float32),
            "layer_ids": np.array(all_layers, dtype=np.int64),
            "labels": np.array(all_labels, dtype=np.float32),
            "track_ids": np.array(all_track_ids, dtype=np.int64),
            "true_tracks": true_tracks,
        })

    return events
