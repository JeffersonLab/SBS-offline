"""
PyTorch Dataset classes for GEM hit denoising and track finding.
"""

import numpy as np
import torch
from torch.utils.data import Dataset

from .root_io import hits_to_graph, generate_synthetic_gem_data, HIT_FEATURE_NAMES


class GEMHitDataset(Dataset):
    """Dataset for autoencoder-based GEM hit denoising.

    Each sample is a fixed-size feature vector for a single detector hit.
    For training the denoiser, clean (signal) hits are used as targets
    and artificially corrupted versions as inputs.

    Parameters
    ----------
    events : list[dict]
        Output of ``generate_synthetic_gem_data`` or equivalent.
        Each dict has "hits" (n_hits, n_feat) and "labels" (n_hits,).
    noise_fraction : float
        Fraction of input features to corrupt during training.
    normalize : bool
        Whether to z-score normalize features.
    signal_only : bool
        If True, only include signal hits (for autoencoder pretraining).
    """

    def __init__(self, events, noise_fraction=0.2, normalize=True, signal_only=False):
        self.noise_fraction = noise_fraction

        all_hits = []
        all_labels = []

        for ev in events:
            hits = ev["hits"]
            labels = ev["labels"]
            if signal_only:
                mask = labels == 1
                hits = hits[mask]
                labels = labels[mask]
            all_hits.append(hits)
            all_labels.append(labels)

        self.hits = np.concatenate(all_hits, axis=0)
        self.labels = np.concatenate(all_labels, axis=0)

        # Compute normalization statistics
        self.normalize = normalize
        if normalize and len(self.hits) > 0:
            self.mean = self.hits.mean(axis=0)
            self.std = self.hits.std(axis=0)
            self.std[self.std < 1e-8] = 1.0
        else:
            self.mean = np.zeros(self.hits.shape[1])
            self.std = np.ones(self.hits.shape[1])

    def __len__(self):
        return len(self.hits)

    def __getitem__(self, idx):
        hit = (self.hits[idx] - self.mean) / self.std
        clean = torch.tensor(hit, dtype=torch.float32)

        # Add corruption noise for denoising autoencoder training
        noisy = clean.clone()
        mask = torch.rand_like(clean) < self.noise_fraction
        noisy[mask] += torch.randn_like(noisy[mask]) * 0.5

        label = torch.tensor(self.labels[idx], dtype=torch.float32)
        return noisy, clean, label

    def get_normalization_params(self):
        return {"mean": self.mean.copy(), "std": self.std.copy()}


class GEMTrackGraphDataset(Dataset):
    """Dataset for GNN-based track finding.

    Each sample is a graph representing one event, where nodes are hits
    and edges connect geometrically compatible hits in adjacent layers.

    Edge labels indicate whether both connected hits belong to the same
    true track (for supervised training).

    Parameters
    ----------
    events : list[dict]
        Output of ``generate_synthetic_gem_data`` or equivalent.
    r_connect : float
        Maximum transverse distance for edge construction (mm).
    k_nearest : int
        Maximum nearest neighbors per hit for edge construction.
    normalize : bool
        Whether to z-score normalize node features.
    """

    def __init__(self, events, r_connect=50.0, k_nearest=8, normalize=True):
        self.graphs = []

        # Compute global normalization from all hits
        all_hits = np.concatenate([ev["hits"] for ev in events], axis=0)
        if normalize and len(all_hits) > 0:
            self.mean = all_hits.mean(axis=0)
            self.std = all_hits.std(axis=0)
            self.std[self.std < 1e-8] = 1.0
        else:
            self.mean = np.zeros(all_hits.shape[1])
            self.std = np.ones(all_hits.shape[1])

        for ev in events:
            hits = ev["hits"]
            layer_ids = ev["layer_ids"]
            track_ids = ev["track_ids"]
            labels = ev["labels"]

            if len(hits) == 0:
                continue

            # Normalize node features
            node_features = (hits - self.mean) / self.std

            # Build graph
            edge_index, edge_attr = hits_to_graph(
                hits, layer_ids, r_connect=r_connect, k_nearest=k_nearest
            )

            # Compute edge truth labels:
            # An edge is "true" if both endpoints are signal hits on the same track
            n_edges = edge_index.shape[1]
            edge_labels = np.zeros(n_edges, dtype=np.float32)
            for e in range(n_edges):
                src, dst = edge_index[0, e], edge_index[1, e]
                if (labels[src] == 1 and labels[dst] == 1 and
                        track_ids[src] >= 0 and track_ids[src] == track_ids[dst]):
                    edge_labels[e] = 1.0

            self.graphs.append({
                "node_features": torch.tensor(node_features, dtype=torch.float32),
                "edge_index": torch.tensor(edge_index, dtype=torch.long),
                "edge_attr": torch.tensor(edge_attr, dtype=torch.float32),
                "edge_labels": torch.tensor(edge_labels, dtype=torch.float32),
                "node_labels": torch.tensor(labels, dtype=torch.float32),
                "track_ids": torch.tensor(track_ids, dtype=torch.long),
                "layer_ids": torch.tensor(layer_ids, dtype=torch.long),
            })

    def __len__(self):
        return len(self.graphs)

    def __getitem__(self, idx):
        return self.graphs[idx]

    def get_normalization_params(self):
        return {"mean": self.mean.copy(), "std": self.std.copy()}


def collate_graphs(batch):
    """Custom collate function for batching variable-size graphs.

    Merges multiple graphs into a single disconnected graph with
    offset node indices, following the PyG convention.
    """
    node_features = []
    edge_index = []
    edge_attr = []
    edge_labels = []
    node_labels = []
    track_ids = []
    layer_ids = []
    batch_idx = []

    node_offset = 0
    for i, graph in enumerate(batch):
        n_nodes = graph["node_features"].shape[0]

        node_features.append(graph["node_features"])
        edge_index.append(graph["edge_index"] + node_offset)
        edge_attr.append(graph["edge_attr"])
        edge_labels.append(graph["edge_labels"])
        node_labels.append(graph["node_labels"])
        track_ids.append(graph["track_ids"])
        layer_ids.append(graph["layer_ids"])
        batch_idx.append(torch.full((n_nodes,), i, dtype=torch.long))

        node_offset += n_nodes

    return {
        "node_features": torch.cat(node_features, dim=0),
        "edge_index": torch.cat(edge_index, dim=1),
        "edge_attr": torch.cat(edge_attr, dim=0),
        "edge_labels": torch.cat(edge_labels, dim=0),
        "node_labels": torch.cat(node_labels, dim=0),
        "track_ids": torch.cat(track_ids, dim=0),
        "layer_ids": torch.cat(layer_ids, dim=0),
        "batch": torch.cat(batch_idx, dim=0),
    }
