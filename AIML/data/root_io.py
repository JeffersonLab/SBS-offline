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
