#!/usr/bin/env python3
"""
Training script for the GEM hit denoising autoencoder.

Usage:
    # Train on synthetic data (for development/testing):
    python train_autoencoder.py --synthetic --epochs 50

    # Train on ROOT file data:
    python train_autoencoder.py --input /path/to/data.root --epochs 100

    # Resume from checkpoint:
    python train_autoencoder.py --input /path/to/data.root --resume checkpoint.pt
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

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.autoencoder import GEMHitAutoencoder, DenoisingLoss
from data.gem_dataset import GEMHitDataset
from data.root_io import generate_synthetic_gem_data


def parse_args():
    parser = argparse.ArgumentParser(description="Train GEM hit denoising autoencoder")
    parser.add_argument("--input", type=str, default=None,
                        help="Path to input ROOT file")
    parser.add_argument("--tree-name", type=str, default="T",
                        help="ROOT tree name")
    parser.add_argument("--branch-prefix", type=str, default="sbs.gems",
                        help="Branch prefix for GEM hits")
    parser.add_argument("--synthetic", action="store_true",
                        help="Use synthetic data for training")
    parser.add_argument("--n-events", type=int, default=5000,
                        help="Number of synthetic events")
    parser.add_argument("--epochs", type=int, default=100,
                        help="Number of training epochs")
    parser.add_argument("--batch-size", type=int, default=512,
                        help="Training batch size")
    parser.add_argument("--lr", type=float, default=1e-3,
                        help="Learning rate")
    parser.add_argument("--latent-dim", type=int, default=8,
                        help="Latent space dimension")
    parser.add_argument("--hidden-dims", type=int, nargs="+", default=[64, 32],
                        help="Hidden layer dimensions")
    parser.add_argument("--noise-fraction", type=float, default=0.2,
                        help="Input corruption fraction")
    parser.add_argument("--recon-weight", type=float, default=1.0,
                        help="Reconstruction loss weight")
    parser.add_argument("--class-weight", type=float, default=0.5,
                        help="Classification loss weight")
    parser.add_argument("--output-dir", type=str, default="checkpoints/autoencoder",
                        help="Output directory for checkpoints")
    parser.add_argument("--resume", type=str, default=None,
                        help="Path to checkpoint to resume from")
    parser.add_argument("--val-fraction", type=float, default=0.15,
                        help="Validation set fraction")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--device", type=str, default=None,
                        help="Device (cuda/cpu, auto-detected if not set)")
    return parser.parse_args()


def train_epoch(model, loader, criterion, optimizer, device):
    model.train()
    total_loss_dict = {}
    n_batches = 0

    for noisy, clean, labels in loader:
        noisy = noisy.to(device)
        clean = clean.to(device)
        labels = labels.to(device)

        optimizer.zero_grad()
        x_recon, z, signal_score = model(noisy)
        loss, loss_dict = criterion(x_recon, clean, signal_score, labels, z)
        loss.backward()

        # Gradient clipping
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        for k, v in loss_dict.items():
            total_loss_dict[k] = total_loss_dict.get(k, 0) + v
        n_batches += 1

    return {k: v / n_batches for k, v in total_loss_dict.items()}


@torch.no_grad()
def validate(model, loader, criterion, device):
    model.eval()
    total_loss_dict = {}
    all_preds = []
    all_labels = []
    n_batches = 0

    for noisy, clean, labels in loader:
        noisy = noisy.to(device)
        clean = clean.to(device)
        labels = labels.to(device)

        x_recon, z, signal_score = model(noisy)
        _, loss_dict = criterion(x_recon, clean, signal_score, labels, z)

        all_preds.append(torch.sigmoid(signal_score).cpu())
        all_labels.append(labels.cpu())

        for k, v in loss_dict.items():
            total_loss_dict[k] = total_loss_dict.get(k, 0) + v
        n_batches += 1

    avg_loss = {k: v / n_batches for k, v in total_loss_dict.items()}

    # Classification metrics
    preds = torch.cat(all_preds)
    labels = torch.cat(all_labels)
    binary_preds = (preds > 0.5).float()
    accuracy = (binary_preds == labels).float().mean().item()

    signal_mask = labels == 1
    noise_mask = labels == 0
    if signal_mask.sum() > 0:
        signal_eff = binary_preds[signal_mask].mean().item()
    else:
        signal_eff = 0.0
    if noise_mask.sum() > 0:
        noise_rej = (1 - binary_preds[noise_mask]).mean().item()
    else:
        noise_rej = 0.0

    avg_loss["accuracy"] = accuracy
    avg_loss["signal_efficiency"] = signal_eff
    avg_loss["noise_rejection"] = noise_rej

    return avg_loss


def main():
    args = parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Device
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
    elif args.input:
        from data.root_io import load_gem_hits_from_root
        print(f"Loading data from {args.input}...")
        hits_list, labels_list = load_gem_hits_from_root(args.input)
        events = []
        for h, l in zip(hits_list, labels_list or [None] * len(hits_list)):
            ev = {"hits": h, "labels": l if l is not None else np.ones(len(h))}
            events.append(ev)
    else:
        print("Error: specify --input or --synthetic")
        sys.exit(1)

    # Create dataset
    dataset = GEMHitDataset(events, noise_fraction=args.noise_fraction)
    print(f"Dataset: {len(dataset)} hits, {dataset.hits.shape[1]} features")

    n_val = int(len(dataset) * args.val_fraction)
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val])

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, drop_last=True)
    val_loader = DataLoader(val_set, batch_size=args.batch_size, shuffle=False,
                            num_workers=0)

    # Model
    n_features = dataset.hits.shape[1]
    model = GEMHitAutoencoder(
        n_features=n_features,
        latent_dim=args.latent_dim,
        hidden_dims=args.hidden_dims,
    ).to(device)

    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Loss and optimizer
    criterion = DenoisingLoss(
        recon_weight=args.recon_weight,
        class_weight=args.class_weight,
    )
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

    # Output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Save normalization parameters
    norm_params = dataset.get_normalization_params()
    np.savez(os.path.join(args.output_dir, "normalization.npz"),
             mean=norm_params["mean"], std=norm_params["std"])

    # Training loop
    history = []
    print(f"\nTraining for {args.epochs} epochs...")
    print("-" * 80)

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
                f"Train loss: {train_loss['total_loss']:.4f} | "
                f"Val loss: {val_loss['total_loss']:.4f} | "
                f"Acc: {val_loss['accuracy']:.3f} | "
                f"Sig.Eff: {val_loss['signal_efficiency']:.3f} | "
                f"Noise.Rej: {val_loss['noise_rejection']:.3f}"
            )

        # Save best model
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

    # Save training history
    with open(os.path.join(args.output_dir, "training_history.json"), "w") as f:
        json.dump(history, f, indent=2)

    print("-" * 80)
    print(f"Training complete. Best val loss: {best_val_loss:.4f}")
    print(f"Models saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
