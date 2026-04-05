"""
Visualization utilities for AIML training and evaluation.
"""

import json
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False


def _check_mpl():
    if not _HAS_MPL:
        raise ImportError("matplotlib is required for visualization. "
                          "Install with: pip install matplotlib")


def plot_event_display(hits, labels=None, track_ids=None, predicted_tracks=None,
                       title="Event Display", output_path=None):
    """Plot a 2D/3D event display of GEM hits.

    Parameters
    ----------
    hits : np.ndarray
        Shape (n_hits, n_features). First 3 columns are x, y, z.
    labels : np.ndarray, optional
        Signal/noise labels for coloring.
    track_ids : np.ndarray, optional
        True track assignments.
    predicted_tracks : list[list[int]], optional
        Predicted track candidates to draw.
    title : str
        Plot title.
    output_path : str, optional
        Save path. If None, shows interactively.
    """
    _check_mpl()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    x, y, z = hits[:, 0], hits[:, 1], hits[:, 2]

    # XZ view
    ax = axes[0]
    if labels is not None:
        colors = np.where(labels == 1, "C0", "C3")
        ax.scatter(z[labels == 0], x[labels == 0], c="C3", s=10, alpha=0.3, label="Noise")
        ax.scatter(z[labels == 1], x[labels == 1], c="C0", s=20, alpha=0.7, label="Signal")
    else:
        ax.scatter(z, x, s=10, alpha=0.5)
    ax.set_xlabel("Z (mm)")
    ax.set_ylabel("X (mm)")
    ax.set_title("X-Z View")
    ax.legend()

    # YZ view
    ax = axes[1]
    if labels is not None:
        ax.scatter(z[labels == 0], y[labels == 0], c="C3", s=10, alpha=0.3, label="Noise")
        ax.scatter(z[labels == 1], y[labels == 1], c="C0", s=20, alpha=0.7, label="Signal")
    else:
        ax.scatter(z, y, s=10, alpha=0.5)
    ax.set_xlabel("Z (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title("Y-Z View")
    ax.legend()

    # XY view
    ax = axes[2]
    if track_ids is not None:
        unique_tracks = np.unique(track_ids[track_ids >= 0])
        cmap = plt.cm.tab10
        for t in unique_tracks:
            mask = track_ids == t
            ax.scatter(x[mask], y[mask], c=[cmap(int(t) % 10)], s=30,
                       alpha=0.8, label=f"Track {t}")
        noise_mask = track_ids < 0
        ax.scatter(x[noise_mask], y[noise_mask], c="gray", s=5, alpha=0.2, label="Noise")
    elif labels is not None:
        ax.scatter(x[labels == 0], y[labels == 0], c="C3", s=10, alpha=0.3, label="Noise")
        ax.scatter(x[labels == 1], y[labels == 1], c="C0", s=20, alpha=0.7, label="Signal")
    else:
        ax.scatter(x, y, s=10, alpha=0.5)

    # Draw predicted tracks
    if predicted_tracks:
        for tidx, track in enumerate(predicted_tracks):
            track_hits = hits[track]
            # Sort by z and draw line
            order = np.argsort(track_hits[:, 2])
            ax.plot(track_hits[order, 0], track_hits[order, 1],
                    "--", linewidth=1.5, alpha=0.7, label=f"Pred {tidx}")

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title("X-Y View")
    ax.legend(fontsize=7, loc="upper right")

    fig.suptitle(title)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.close(fig)
    return fig


def plot_training_curves(history_path, output_path=None):
    """Plot training and validation loss curves from a training history JSON.

    Parameters
    ----------
    history_path : str
        Path to training_history.json.
    output_path : str, optional
        Save path.
    """
    _check_mpl()

    with open(history_path) as f:
        history = json.load(f)

    epochs = [h["epoch"] for h in history]
    train_loss = [h["train"]["total_loss"] for h in history]
    val_loss = [h["val"]["total_loss"] for h in history]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.plot(epochs, train_loss, label="Train")
    ax.plot(epochs, val_loss, label="Validation")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Total Loss")
    ax.set_title("Training Curves")
    ax.legend()
    ax.set_yscale("log")

    # Check for classification metrics
    ax = axes[1]
    keys_to_plot = []
    for key in ["accuracy", "signal_efficiency", "noise_rejection",
                "edge_f1", "edge_precision", "edge_recall", "node_accuracy"]:
        if key in history[0]["val"]:
            keys_to_plot.append(key)

    for key in keys_to_plot:
        values = [h["val"][key] for h in history]
        ax.plot(epochs, values, label=key)

    ax.set_xlabel("Epoch")
    ax.set_ylabel("Metric")
    ax.set_title("Validation Metrics")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.close(fig)
    return fig


def plot_roc_curve(predictions, labels, title="ROC Curve", output_path=None):
    """Plot ROC curve for binary classification.

    Parameters
    ----------
    predictions : np.ndarray
        Predicted probabilities.
    labels : np.ndarray
        True binary labels.
    """
    _check_mpl()

    thresholds = np.linspace(0, 1, 200)
    tpr_list = []
    fpr_list = []

    for t in thresholds:
        pred_pos = predictions >= t
        tp = ((pred_pos) & (labels == 1)).sum()
        fp = ((pred_pos) & (labels == 0)).sum()
        fn = ((~pred_pos) & (labels == 1)).sum()
        tn = ((~pred_pos) & (labels == 0)).sum()

        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
        tpr_list.append(tpr)
        fpr_list.append(fpr)

    # AUC via trapezoidal rule
    fpr_arr = np.array(fpr_list)
    tpr_arr = np.array(tpr_list)
    sorted_idx = np.argsort(fpr_arr)
    auc = np.trapz(tpr_arr[sorted_idx], fpr_arr[sorted_idx])

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(fpr_arr, tpr_arr, linewidth=2, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="gray", linewidth=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.close(fig)
    return fig
