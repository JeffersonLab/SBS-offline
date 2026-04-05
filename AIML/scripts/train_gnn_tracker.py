#!/usr/bin/env python3
"""
Training script for the GNN-based track finder.

Usage:
    # Train on synthetic data:
    python train_gnn_tracker.py --synthetic --epochs 50

    # Train on ROOT file:
    python train_gnn_tracker.py --input /path/to/data.root --epochs 100
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.gnn_tracker import TrackGNN, TrackGNNLoss, extract_tracks_from_edges
from data.gem_dataset import GEMTrackGraphDataset, collate_graphs
from data.root_io import generate_synthetic_gem_data


def parse_args():
    parser = argparse.ArgumentParser(description="Train GNN track finder")
    parser.add_argument("--input", type=str, default=None,
                        help="Input ROOT file path")
    parser.add_argument("--synthetic", action="store_true",
                        help="Use synthetic data")
    parser.add_argument("--n-events", type=int, default=2000,
                        help="Number of synthetic events")
    parser.add_argument("--epochs", type=int, default=80,
                        help="Training epochs")
    parser.add_argument("--batch-size", type=int, default=16,
                        help="Batch size (number of events/graphs)")
    parser.add_argument("--lr", type=float, default=5e-4,
                        help="Learning rate")
    parser.add_argument("--hidden-dim", type=int, default=64,
                        help="Hidden dimension")
    parser.add_argument("--n-iterations", type=int, default=4,
                        help="Number of message-passing iterations")
    parser.add_argument("--r-connect", type=float, default=50.0,
                        help="Edge connection radius (mm)")
    parser.add_argument("--k-nearest", type=int, default=8,
                        help="Max nearest neighbors per hit")
    parser.add_argument("--pos-weight", type=float, default=5.0,
                        help="Positive edge weight for BCE loss")
    parser.add_argument("--output-dir", type=str, default="checkpoints/gnn_tracker",
                        help="Output directory")
    parser.add_argument("--resume", type=str, default=None,
                        help="Resume from checkpoint")
    parser.add_argument("--val-fraction", type=float, default=0.15,
                        help="Validation fraction")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", type=str, default=None)
    return parser.parse_args()


def train_epoch(model, loader, criterion, optimizer, device):
    model.train()
    total_loss_dict = {}
    n_batches = 0

    for batch in loader:
        node_feat = batch["node_features"].to(device)
        edge_idx = batch["edge_index"].to(device)
        edge_attr = batch["edge_attr"].to(device)
        edge_labels = batch["edge_labels"].to(device)
        node_labels = batch["node_labels"].to(device)

        if edge_idx.shape[1] == 0:
            continue

        optimizer.zero_grad()
        edge_scores, node_scores, inter_scores = model(node_feat, edge_idx, edge_attr)
        loss, loss_dict = criterion(
            edge_scores, node_scores, inter_scores, edge_labels, node_labels
        )
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        for k, v in loss_dict.items():
            total_loss_dict[k] = total_loss_dict.get(k, 0) + v
        n_batches += 1

    if n_batches == 0:
        return {"total_loss": 0.0}
    return {k: v / n_batches for k, v in total_loss_dict.items()}


@torch.no_grad()
def validate(model, loader, criterion, device):
    model.eval()
    total_loss_dict = {}
    all_edge_preds = []
    all_edge_labels = []
    all_node_preds = []
    all_node_labels = []
    n_batches = 0

    for batch in loader:
        node_feat = batch["node_features"].to(device)
        edge_idx = batch["edge_index"].to(device)
        edge_attr = batch["edge_attr"].to(device)
        edge_labels = batch["edge_labels"].to(device)
        node_labels = batch["node_labels"].to(device)

        if edge_idx.shape[1] == 0:
            continue

        edge_scores, node_scores, inter_scores = model(node_feat, edge_idx, edge_attr)
        _, loss_dict = criterion(
            edge_scores, node_scores, inter_scores, edge_labels, node_labels
        )

        all_edge_preds.append(torch.sigmoid(edge_scores).cpu())
        all_edge_labels.append(edge_labels.cpu())
        all_node_preds.append(torch.sigmoid(node_scores).cpu())
        all_node_labels.append(node_labels.cpu())

        for k, v in loss_dict.items():
            total_loss_dict[k] = total_loss_dict.get(k, 0) + v
        n_batches += 1

    if n_batches == 0:
        return {"total_loss": 0.0}

    avg_loss = {k: v / n_batches for k, v in total_loss_dict.items()}

    # Edge metrics
    edge_preds = torch.cat(all_edge_preds)
    edge_labels = torch.cat(all_edge_labels)
    edge_binary = (edge_preds > 0.5).float()

    true_pos = ((edge_binary == 1) & (edge_labels == 1)).sum().item()
    false_pos = ((edge_binary == 1) & (edge_labels == 0)).sum().item()
    false_neg = ((edge_binary == 0) & (edge_labels == 1)).sum().item()

    precision = true_pos / (true_pos + false_pos) if (true_pos + false_pos) > 0 else 0
    recall = true_pos / (true_pos + false_neg) if (true_pos + false_neg) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_loss["edge_precision"] = precision
    avg_loss["edge_recall"] = recall
    avg_loss["edge_f1"] = f1

    # Node metrics
    node_preds = torch.cat(all_node_preds)
    node_labels = torch.cat(all_node_labels)
    node_acc = ((node_preds > 0.5).float() == node_labels).float().mean().item()
    avg_loss["node_accuracy"] = node_acc

    return avg_loss


def evaluate_tracking(model, dataset, device, threshold=0.5):
    """Evaluate full tracking pipeline on individual events."""
    model.eval()
    efficiencies = []
    fake_rates = []

    for i in range(min(len(dataset), 100)):
        graph = dataset[i]
        node_feat = graph["node_features"].unsqueeze(0).to(device) if graph["node_features"].dim() == 2 else graph["node_features"].to(device)
        edge_idx = graph["edge_index"].to(device)
        edge_attr = graph["edge_attr"].to(device)
        track_ids = graph["track_ids"]

        if edge_idx.shape[1] == 0:
            continue

        with torch.no_grad():
            edge_scores, _, _ = model(
                graph["node_features"].to(device), edge_idx, edge_attr
            )
            edge_probs = torch.sigmoid(edge_scores)

        # Extract track candidates
        candidates = extract_tracks_from_edges(edge_idx, edge_probs, threshold)

        # Evaluate: how many true tracks are found?
        true_track_ids = set(track_ids[track_ids >= 0].numpy().tolist())
        n_true_tracks = len(true_track_ids)

        found_tracks = 0
        fake_tracks = 0
        for candidate in candidates:
            cand_track_ids = track_ids[candidate]
            valid = cand_track_ids[cand_track_ids >= 0]
            if len(valid) == 0:
                fake_tracks += 1
                continue
            # Majority vote
            majority_id = np.bincount(valid.numpy()).argmax()
            purity = (valid == majority_id).float().mean().item()
            if purity > 0.5 and len(valid) >= 3:
                found_tracks += 1
            else:
                fake_tracks += 1

        eff = found_tracks / n_true_tracks if n_true_tracks > 0 else 0
        fr = fake_tracks / len(candidates) if candidates else 0
        efficiencies.append(eff)
        fake_rates.append(fr)

    return {
        "tracking_efficiency": np.mean(efficiencies) if efficiencies else 0,
        "fake_rate": np.mean(fake_rates) if fake_rates else 0,
    }


def main():
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load data
    if args.synthetic:
        print(f"Generating {args.n_events} synthetic events...")
        events = generate_synthetic_gem_data(
            n_events=args.n_events, n_layers=5,
            n_signal_tracks=2, n_noise_hits_per_layer=10,
        )
    else:
        print("Error: ROOT input for GNN not yet implemented. Use --synthetic.")
        sys.exit(1)

    # Build graph dataset
    print("Building graph dataset...")
    dataset = GEMTrackGraphDataset(
        events, r_connect=args.r_connect, k_nearest=args.k_nearest
    )
    print(f"Dataset: {len(dataset)} event graphs")

    n_val = int(len(dataset) * args.val_fraction)
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val])

    train_loader = DataLoader(
        train_set, batch_size=args.batch_size, shuffle=True,
        collate_fn=collate_graphs, num_workers=0, drop_last=True,
    )
    val_loader = DataLoader(
        val_set, batch_size=args.batch_size, shuffle=False,
        collate_fn=collate_graphs, num_workers=0,
    )

    # Model
    sample = dataset[0]
    node_dim = sample["node_features"].shape[1]
    edge_dim = sample["edge_attr"].shape[1]

    model = TrackGNN(
        node_dim=node_dim, edge_dim=edge_dim,
        hidden_dim=args.hidden_dim, n_iterations=args.n_iterations,
    ).to(device)

    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Loss and optimizer
    criterion = TrackGNNLoss(pos_weight=args.pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    # Resume
    start_epoch = 0
    best_val_loss = float("inf")
    if args.resume:
        ckpt = torch.load(args.resume, map_location=device)
        model.load_state_dict(ckpt["model_state_dict"])
        optimizer.load_state_dict(ckpt["optimizer_state_dict"])
        start_epoch = ckpt.get("epoch", 0) + 1
        best_val_loss = ckpt.get("best_val_loss", float("inf"))
        print(f"Resumed from epoch {start_epoch}")

    os.makedirs(args.output_dir, exist_ok=True)

    # Save normalization params
    norm_params = dataset.get_normalization_params()
    np.savez(os.path.join(args.output_dir, "normalization.npz"),
             mean=norm_params["mean"], std=norm_params["std"])

    # Training loop
    history = []
    print(f"\nTraining for {args.epochs} epochs...")
    print("-" * 90)

    for epoch in range(start_epoch, args.epochs):
        t0 = time.time()

        train_loss = train_epoch(model, train_loader, criterion, optimizer, device)
        val_loss = validate(model, val_loader, criterion, device)
        scheduler.step()
        dt = time.time() - t0

        history.append({"epoch": epoch, "train": train_loss, "val": val_loss})

        if epoch % 5 == 0 or epoch == args.epochs - 1:
            print(
                f"Epoch {epoch:3d}/{args.epochs} ({dt:.1f}s) | "
                f"Train: {train_loss['total_loss']:.4f} | "
                f"Val: {val_loss['total_loss']:.4f} | "
                f"Edge P/R/F1: {val_loss.get('edge_precision', 0):.3f}/"
                f"{val_loss.get('edge_recall', 0):.3f}/"
                f"{val_loss.get('edge_f1', 0):.3f} | "
                f"Node Acc: {val_loss.get('node_accuracy', 0):.3f}"
            )

        if val_loss["total_loss"] < best_val_loss:
            best_val_loss = val_loss["total_loss"]
            torch.save({
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "best_val_loss": best_val_loss,
                "args": vars(args),
            }, os.path.join(args.output_dir, "best_model.pt"))

    # Save final model
    torch.save({
        "epoch": args.epochs - 1,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "best_val_loss": best_val_loss,
        "args": vars(args),
    }, os.path.join(args.output_dir, "final_model.pt"))

    # Full tracking evaluation
    print("\nEvaluating tracking performance...")
    tracking_metrics = evaluate_tracking(model, dataset, device)
    print(f"Tracking efficiency: {tracking_metrics['tracking_efficiency']:.3f}")
    print(f"Fake rate: {tracking_metrics['fake_rate']:.3f}")

    # Save history
    with open(os.path.join(args.output_dir, "training_history.json"), "w") as f:
        json.dump(history, f, indent=2)

    print("-" * 90)
    print(f"Training complete. Best val loss: {best_val_loss:.4f}")
    print(f"Models saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
