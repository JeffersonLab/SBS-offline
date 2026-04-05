# SBS-offline AIML Module

Machine learning components for the SuperBigBite (SBS) detector reconstruction at Jefferson Lab. This module provides two main capabilities:

1. **Autoencoder-based hit denoising** - Removes electronic noise and background hits from GEM detector data using a denoising autoencoder with joint signal/noise classification.

2. **GNN-based track finding** - Finds particle tracks through GEM detector layers using a Graph Neural Network that classifies hit-pair connections.

## Architecture

### Autoencoder Denoiser (`models/autoencoder.py`)

- **GEMHitAutoencoder**: Fully-connected denoising autoencoder operating on per-hit feature vectors (position, energy, timing, correlations). Maps through a compressed latent space that filters noise while preserving signal. Includes a classifier head for signal/noise discrimination.

- **GEMConvAutoencoder**: 1D convolutional autoencoder for raw ADC strip waveforms (6 time samples). Denoises waveform shapes before clustering.

### GNN Track Finder (`models/gnn_tracker.py`)

- **TrackGNN**: Interaction-network-style GNN with iterative message passing. Hits are nodes, edges connect geometrically compatible hits in adjacent layers. After N message-passing iterations with edge classification and node updates, final edge scores indicate same-track probabilities.

- **EdgeClassifier**: Lightweight MLP-based edge classifier for fast initial filtering before full GNN inference.

### Data Pipeline (`data/`)

- **root_io.py**: Reads GEM hit/track data from SBS-offline ROOT files via `uproot`. Also provides synthetic data generation for development.
- **gem_dataset.py**: PyTorch Dataset classes for both autoencoder training (per-hit) and GNN training (per-event graphs).

## Quick Start

### Install dependencies

```bash
pip install -r requirements.txt
```

### Train on synthetic data

```bash
# Train autoencoder
python scripts/train_autoencoder.py --synthetic --epochs 50 --output-dir checkpoints/autoencoder

# Train GNN tracker
python scripts/train_gnn_tracker.py --synthetic --epochs 50 --output-dir checkpoints/gnn_tracker
```

### Run inference

```bash
python scripts/inference.py \
    --synthetic \
    --autoencoder-checkpoint checkpoints/autoencoder/best_model.pt \
    --gnn-checkpoint checkpoints/gnn_tracker/best_model.pt \
    --output results.npz
```

### Train on real data

```bash
# From SBS-offline replay ROOT files:
python scripts/train_autoencoder.py \
    --input /path/to/replay_output.root \
    --branch-prefix sbs.gems \
    --epochs 100
```

## Integration with SBS-offline

The AIML module reads from and writes to the same ROOT tree structure used by the existing SBS-offline reconstruction. Key mappings:

| SBS-offline Structure | AIML Feature |
|---|---|
| `sbsgemhit_t.xghit/yghit/zghit` | Node position (x, y, z) |
| `sbsgemhit_t.Ehit` | Hit energy |
| `sbsgemhit_t.thit/thitcorr` | Hit timing |
| `sbsgemhit_t.ADCasym/tdiff` | U/V correlation features |
| `sbsgemhit_t.corrcoeff_clust/strip` | Cluster quality |
| `sbsgemhit_t.ontrack` | Training label (MC truth) |
| `SBSGEMTrackerBase` track parameters | Validation reference |

## Directory Structure

```
AIML/
├── __init__.py
├── README.md
├── requirements.txt
├── models/
│   ├── __init__.py
│   ├── autoencoder.py      # Denoising autoencoder models
│   └── gnn_tracker.py      # GNN track finding models
├── data/
│   ├── __init__.py
│   ├── root_io.py          # ROOT file I/O and synthetic data
│   └── gem_dataset.py      # PyTorch datasets
├── utils/
│   ├── __init__.py
│   ├── metrics.py          # Tracking efficiency, fake rate metrics
│   └── visualization.py    # Event displays, training curves
└── scripts/
    ├── __init__.py
    ├── train_autoencoder.py   # Autoencoder training script
    ├── train_gnn_tracker.py   # GNN training script
    └── inference.py           # Full inference pipeline
```
