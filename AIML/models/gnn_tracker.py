"""
Graph Neural Network models for track finding in the SBS GEM detector.

Uses a message-passing architecture where:
- Nodes represent detector hits
- Edges connect hits in adjacent detector layers
- Edge classification determines which hit pairs belong to the same track

The approach follows the Interaction Network / Edge Network paradigm
used successfully in HEP tracking (TrackML, ATLAS, CMS).

This implementation uses only PyTorch (no PyG dependency required),
making it easier to deploy in the JLab computing environment.
"""

import torch
import torch.nn as nn


class EdgeNetwork(nn.Module):
    """Classifies edges as true (same-track) or fake based on
    node features of the connected hits and edge geometric features.

    Parameters
    ----------
    node_dim : int
        Dimension of node feature vectors.
    edge_dim : int
        Dimension of edge feature vectors.
    hidden_dim : int
        Hidden layer size.
    """

    def __init__(self, node_dim, edge_dim, hidden_dim=64):
        super().__init__()
        # Input: concatenation of [src_features, dst_features, edge_features]
        input_dim = 2 * node_dim + edge_dim
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, node_features, edge_index, edge_attr):
        """Compute edge scores.

        Parameters
        ----------
        node_features : torch.Tensor
            Shape (n_nodes, node_dim).
        edge_index : torch.Tensor
            Shape (2, n_edges).
        edge_attr : torch.Tensor
            Shape (n_edges, edge_dim).

        Returns
        -------
        edge_scores : torch.Tensor
            Shape (n_edges,) - logits for edge classification.
        """
        src = node_features[edge_index[0]]
        dst = node_features[edge_index[1]]
        edge_input = torch.cat([src, dst, edge_attr], dim=-1)
        return self.net(edge_input).squeeze(-1)


class NodeNetwork(nn.Module):
    """Updates node features by aggregating messages from neighboring nodes,
    weighted by edge scores.

    Parameters
    ----------
    node_dim : int
        Dimension of node feature vectors.
    edge_dim : int
        Dimension of edge feature vectors.
    hidden_dim : int
        Hidden layer size.
    """

    def __init__(self, node_dim, edge_dim, hidden_dim=64):
        super().__init__()
        # Message MLP: transforms [src, edge_attr] into a message
        self.message_mlp = nn.Sequential(
            nn.Linear(node_dim + edge_dim, hidden_dim),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_dim, node_dim),
        )
        # Update MLP: combines old node features with aggregated messages
        self.update_mlp = nn.Sequential(
            nn.Linear(2 * node_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_dim, node_dim),
        )

    def forward(self, node_features, edge_index, edge_attr, edge_weights):
        """Update node features via message passing.

        Parameters
        ----------
        node_features : torch.Tensor
            Shape (n_nodes, node_dim).
        edge_index : torch.Tensor
            Shape (2, n_edges).
        edge_attr : torch.Tensor
            Shape (n_edges, edge_dim).
        edge_weights : torch.Tensor
            Shape (n_edges,) - sigmoid edge scores used as attention weights.

        Returns
        -------
        updated_features : torch.Tensor
            Shape (n_nodes, node_dim).
        """
        src_idx, dst_idx = edge_index[0], edge_index[1]
        src_features = node_features[src_idx]

        # Compute messages
        msg_input = torch.cat([src_features, edge_attr], dim=-1)
        messages = self.message_mlp(msg_input)

        # Weight messages by edge confidence
        messages = messages * edge_weights.unsqueeze(-1)

        # Aggregate messages at destination nodes (sum)
        n_nodes = node_features.shape[0]
        aggregated = torch.zeros_like(node_features)
        aggregated.index_add_(0, dst_idx, messages)

        # Update node features
        combined = torch.cat([node_features, aggregated], dim=-1)
        updated = self.update_mlp(combined)

        return updated


class TrackGNN(nn.Module):
    """Full GNN model for track finding via iterative edge classification.

    Alternates between:
    1. Edge classification (which hit pairs are on the same track?)
    2. Node update (refine hit representations using neighbor information)

    After several message-passing iterations, the final edge scores
    indicate track segment probabilities.

    Parameters
    ----------
    node_dim : int
        Input node feature dimension (default: 14 for GEM hits).
    edge_dim : int
        Input edge feature dimension (default: 5 for dx,dy,dz,dr,dphi).
    hidden_dim : int
        Hidden layer dimension.
    n_iterations : int
        Number of message-passing iterations.
    """

    def __init__(self, node_dim=14, edge_dim=5, hidden_dim=64, n_iterations=4):
        super().__init__()

        self.n_iterations = n_iterations

        # Initial node feature transform
        self.node_encoder = nn.Sequential(
            nn.Linear(node_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.LeakyReLU(0.1),
        )

        # Initial edge feature transform
        self.edge_encoder = nn.Sequential(
            nn.Linear(edge_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.LeakyReLU(0.1),
        )

        # Iterative edge and node networks (shared across iterations)
        self.edge_networks = nn.ModuleList([
            EdgeNetwork(hidden_dim, hidden_dim, hidden_dim)
            for _ in range(n_iterations)
        ])
        self.node_networks = nn.ModuleList([
            NodeNetwork(hidden_dim, hidden_dim, hidden_dim)
            for _ in range(n_iterations)
        ])

        # Final edge classifier
        self.edge_classifier = EdgeNetwork(hidden_dim, hidden_dim, hidden_dim)

        # Node classifier (signal vs noise)
        self.node_classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LeakyReLU(0.1),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, node_features, edge_index, edge_attr):
        """Forward pass.

        Parameters
        ----------
        node_features : torch.Tensor
            Shape (n_nodes, node_dim).
        edge_index : torch.Tensor
            Shape (2, n_edges).
        edge_attr : torch.Tensor
            Shape (n_edges, edge_dim).

        Returns
        -------
        edge_scores : torch.Tensor
            Final edge classification logits, shape (n_edges,).
        node_scores : torch.Tensor
            Node signal/noise logits, shape (n_nodes,).
        intermediate_edge_scores : list[torch.Tensor]
            Edge scores at each iteration (for auxiliary loss).
        """
        # Encode initial features
        x = self.node_encoder(node_features)
        e = self.edge_encoder(edge_attr)

        intermediate_edge_scores = []

        # Iterative message passing
        for i in range(self.n_iterations):
            # Edge classification
            edge_scores_i = self.edge_networks[i](x, edge_index, e)
            edge_weights_i = torch.sigmoid(edge_scores_i)
            intermediate_edge_scores.append(edge_scores_i)

            # Node update using edge-weighted messages
            x = x + self.node_networks[i](x, edge_index, e, edge_weights_i)

        # Final edge classification
        final_edge_scores = self.edge_classifier(x, edge_index, e)

        # Node classification
        node_scores = self.node_classifier(x).squeeze(-1)

        return final_edge_scores, node_scores, intermediate_edge_scores


class EdgeClassifier(nn.Module):
    """Lightweight edge classifier for fast inference.

    A simpler model that classifies edges based on concatenated node
    features without iterative message passing. Useful for initial
    edge filtering before running the full GNN.

    Parameters
    ----------
    node_dim : int
        Node feature dimension.
    edge_dim : int
        Edge feature dimension.
    hidden_dim : int
        Hidden layer size.
    n_layers : int
        Number of hidden layers.
    """

    def __init__(self, node_dim=14, edge_dim=5, hidden_dim=64, n_layers=3):
        super().__init__()

        input_dim = 2 * node_dim + edge_dim
        layers = []
        in_dim = input_dim
        for _ in range(n_layers):
            layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.LeakyReLU(0.1),
            ])
            in_dim = hidden_dim
        layers.append(nn.Linear(in_dim, 1))

        self.net = nn.Sequential(*layers)

    def forward(self, node_features, edge_index, edge_attr):
        """Classify edges.

        Returns
        -------
        edge_scores : torch.Tensor
            Edge logits, shape (n_edges,).
        """
        src = node_features[edge_index[0]]
        dst = node_features[edge_index[1]]
        x = torch.cat([src, dst, edge_attr], dim=-1)
        return self.net(x).squeeze(-1)


class TrackGNNLoss(nn.Module):
    """Combined loss for GNN track finding training.

    Includes:
    - Edge classification BCE loss (primary objective)
    - Node classification BCE loss (auxiliary signal/noise)
    - Intermediate iteration losses (deep supervision)

    Parameters
    ----------
    edge_weight : float
        Weight for final edge loss.
    node_weight : float
        Weight for node classification loss.
    intermediate_weight : float
        Weight for intermediate edge losses.
    pos_weight : float
        Positive class weight for edge BCE (handles class imbalance,
        since true edges are typically much rarer than fake edges).
    """

    def __init__(self, edge_weight=1.0, node_weight=0.3,
                 intermediate_weight=0.1, pos_weight=5.0):
        super().__init__()
        self.edge_weight = edge_weight
        self.node_weight = node_weight
        self.intermediate_weight = intermediate_weight
        self.edge_bce = nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor([pos_weight])
        )
        self.node_bce = nn.BCEWithLogitsLoss()

    def forward(self, edge_scores, node_scores, intermediate_scores,
                edge_labels, node_labels):
        """Compute combined loss.

        Returns
        -------
        loss : torch.Tensor
        loss_dict : dict
        """
        # Move pos_weight to same device as scores
        if self.edge_bce.pos_weight.device != edge_scores.device:
            self.edge_bce.pos_weight = self.edge_bce.pos_weight.to(edge_scores.device)

        # Final edge loss
        edge_loss = self.edge_bce(edge_scores, edge_labels)

        # Node loss
        node_loss = self.node_bce(node_scores, node_labels)

        # Intermediate edge losses (deep supervision)
        inter_loss = torch.tensor(0.0, device=edge_scores.device)
        for scores_i in intermediate_scores:
            inter_loss = inter_loss + self.edge_bce(scores_i, edge_labels)
        if intermediate_scores:
            inter_loss = inter_loss / len(intermediate_scores)

        total = (self.edge_weight * edge_loss +
                 self.node_weight * node_loss +
                 self.intermediate_weight * inter_loss)

        loss_dict = {
            "edge_loss": edge_loss.item(),
            "node_loss": node_loss.item(),
            "intermediate_loss": inter_loss.item(),
            "total_loss": total.item(),
        }

        return total, loss_dict


def extract_tracks_from_edges(edge_index, edge_scores, threshold=0.5):
    """Post-processing: extract track candidates from classified edges.

    Uses connected components on the graph of high-score edges to
    group hits into track candidates.

    Parameters
    ----------
    edge_index : torch.Tensor
        Shape (2, n_edges).
    edge_scores : torch.Tensor
        Edge probabilities (after sigmoid), shape (n_edges,).
    threshold : float
        Minimum edge score to consider an edge as "true".

    Returns
    -------
    track_candidates : list[list[int]]
        List of tracks, each being a list of node indices.
    """
    # Filter edges by threshold
    mask = edge_scores > threshold
    if mask.sum() == 0:
        return []

    filtered_edges = edge_index[:, mask].cpu().numpy()

    # Union-Find for connected components
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(filtered_edges.shape[1]):
        union(int(filtered_edges[0, i]), int(filtered_edges[1, i]))

    # Group nodes by component
    components = {}
    for node in parent:
        root = find(node)
        components.setdefault(root, []).append(node)

    # Also include isolated nodes that appear in edges but aren't in parent
    all_nodes = set(filtered_edges.flatten())
    for node in all_nodes:
        root = find(int(node))
        if root not in components:
            components[root] = [int(node)]
        elif int(node) not in components[root]:
            components[root].append(int(node))

    # Filter out single-hit "tracks"
    track_candidates = [
        sorted(nodes) for nodes in components.values() if len(nodes) >= 3
    ]

    return track_candidates
