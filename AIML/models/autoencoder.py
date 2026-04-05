"""
Autoencoder models for GEM hit denoising in the SBS detector.

Two architectures are provided:

1. GEMHitAutoencoder - Fully-connected denoising autoencoder for per-hit
   feature vectors. Learns a compressed latent representation that filters
   noise while preserving signal hit characteristics.

2. GEMConvAutoencoder - 1D convolutional autoencoder for strip-level ADC
   time-sample waveforms. Denoises the raw ADC samples before clustering.
"""

import torch
import torch.nn as nn


class GEMHitAutoencoder(nn.Module):
    """Denoising autoencoder for GEM 2D hit feature vectors.

    Takes corrupted/noisy hit feature vectors and reconstructs clean versions.
    The bottleneck latent space captures the essential signal characteristics
    while filtering electronic noise, pedestal fluctuations, and background
    hits.

    Architecture:
        Input (n_features) -> Encoder -> Latent (latent_dim) -> Decoder -> Output (n_features)

    The encoder also produces a signal/noise classification score, enabling
    joint denoising + hit classification.

    Parameters
    ----------
    n_features : int
        Number of input hit features (default: 14, matching HIT_FEATURE_NAMES).
    latent_dim : int
        Dimension of the bottleneck latent space.
    hidden_dims : list[int]
        Hidden layer sizes for encoder/decoder.
    dropout : float
        Dropout rate for regularization.
    """

    def __init__(self, n_features=14, latent_dim=8, hidden_dims=None, dropout=0.1):
        super().__init__()

        if hidden_dims is None:
            hidden_dims = [64, 32]

        self.n_features = n_features
        self.latent_dim = latent_dim

        # Encoder
        encoder_layers = []
        in_dim = n_features
        for h_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(in_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.LeakyReLU(0.1),
                nn.Dropout(dropout),
            ])
            in_dim = h_dim
        encoder_layers.append(nn.Linear(in_dim, latent_dim))
        self.encoder = nn.Sequential(*encoder_layers)

        # Decoder (mirror of encoder)
        decoder_layers = []
        in_dim = latent_dim
        for h_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(in_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.LeakyReLU(0.1),
                nn.Dropout(dropout),
            ])
            in_dim = h_dim
        decoder_layers.append(nn.Linear(in_dim, n_features))
        self.decoder = nn.Sequential(*decoder_layers)

        # Signal/noise classifier head (from latent space)
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, 16),
            nn.LeakyReLU(0.1),
            nn.Linear(16, 1),
        )

    def encode(self, x):
        """Encode input to latent representation."""
        return self.encoder(x)

    def decode(self, z):
        """Decode latent representation to reconstructed features."""
        return self.decoder(z)

    def classify(self, z):
        """Classify latent representation as signal (1) or noise (0)."""
        return self.classifier(z).squeeze(-1)

    def forward(self, x):
        """Full forward pass: encode, decode, and classify.

        Parameters
        ----------
        x : torch.Tensor
            Input hit features, shape (batch, n_features).

        Returns
        -------
        x_recon : torch.Tensor
            Reconstructed (denoised) features, shape (batch, n_features).
        z : torch.Tensor
            Latent representation, shape (batch, latent_dim).
        signal_score : torch.Tensor
            Signal probability logit, shape (batch,).
        """
        z = self.encode(x)
        x_recon = self.decode(z)
        signal_score = self.classify(z)
        return x_recon, z, signal_score


class GEMConvAutoencoder(nn.Module):
    """1D convolutional denoising autoencoder for GEM strip ADC waveforms.

    Operates on raw ADC time samples (typically 6 samples per strip)
    across strips in a cluster. Learns to denoise the waveform shape
    before peak-finding and clustering.

    Parameters
    ----------
    n_time_samples : int
        Number of ADC time samples per strip (default: 6 for SBS GEM).
    n_channels : int
        Number of input channels (1 for single-strip, or n_strips for cluster).
    latent_channels : int
        Number of channels in the bottleneck.
    """

    def __init__(self, n_time_samples=6, n_channels=1, latent_channels=4):
        super().__init__()

        self.n_time_samples = n_time_samples

        # Encoder: Conv1D layers to compress temporal information
        self.encoder = nn.Sequential(
            nn.Conv1d(n_channels, 16, kernel_size=3, padding=1),
            nn.BatchNorm1d(16),
            nn.LeakyReLU(0.1),
            nn.Conv1d(16, 32, kernel_size=3, padding=1),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(0.1),
            nn.Conv1d(32, latent_channels, kernel_size=3, padding=1),
            nn.BatchNorm1d(latent_channels),
            nn.LeakyReLU(0.1),
        )

        # Decoder: transpose conv to reconstruct waveform
        self.decoder = nn.Sequential(
            nn.Conv1d(latent_channels, 32, kernel_size=3, padding=1),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(0.1),
            nn.Conv1d(32, 16, kernel_size=3, padding=1),
            nn.BatchNorm1d(16),
            nn.LeakyReLU(0.1),
            nn.Conv1d(16, n_channels, kernel_size=3, padding=1),
        )

    def forward(self, x):
        """Forward pass.

        Parameters
        ----------
        x : torch.Tensor
            Input waveform, shape (batch, n_channels, n_time_samples).

        Returns
        -------
        x_recon : torch.Tensor
            Denoised waveform, same shape as input.
        z : torch.Tensor
            Latent representation, shape (batch, latent_channels, n_time_samples).
        """
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return x_recon, z


class DenoisingLoss(nn.Module):
    """Combined loss for denoising autoencoder training.

    Combines:
    - MSE reconstruction loss (denoising quality)
    - BCE classification loss (signal/noise discrimination)
    - KL divergence regularization on latent space (optional)

    Parameters
    ----------
    recon_weight : float
        Weight for reconstruction loss.
    class_weight : float
        Weight for classification loss.
    kl_weight : float
        Weight for KL divergence regularization.
    """

    def __init__(self, recon_weight=1.0, class_weight=0.5, kl_weight=0.01):
        super().__init__()
        self.recon_weight = recon_weight
        self.class_weight = class_weight
        self.kl_weight = kl_weight
        self.mse = nn.MSELoss()
        self.bce = nn.BCEWithLogitsLoss()

    def forward(self, x_recon, x_target, signal_score, signal_label, z=None):
        """Compute combined loss.

        Parameters
        ----------
        x_recon : torch.Tensor
            Reconstructed features.
        x_target : torch.Tensor
            Clean target features.
        signal_score : torch.Tensor
            Signal classification logits.
        signal_label : torch.Tensor
            True signal/noise labels (1/0).
        z : torch.Tensor, optional
            Latent representation for KL regularization.

        Returns
        -------
        loss : torch.Tensor
            Total combined loss.
        loss_dict : dict
            Individual loss components.
        """
        recon_loss = self.mse(x_recon, x_target)
        class_loss = self.bce(signal_score, signal_label)

        loss = self.recon_weight * recon_loss + self.class_weight * class_loss

        loss_dict = {
            "recon_loss": recon_loss.item(),
            "class_loss": class_loss.item(),
        }

        if z is not None and self.kl_weight > 0:
            # Approximate KL divergence assuming unit Gaussian prior
            kl_loss = 0.5 * torch.mean(z.pow(2))
            loss = loss + self.kl_weight * kl_loss
            loss_dict["kl_loss"] = kl_loss.item()

        loss_dict["total_loss"] = loss.item()
        return loss, loss_dict
