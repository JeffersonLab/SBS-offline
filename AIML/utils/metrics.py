"""
Evaluation metrics for denoising and tracking performance.
"""

import numpy as np


def denoising_metrics(predictions, targets, labels, threshold=0.5):
    """Compute denoising autoencoder performance metrics.

    Parameters
    ----------
    predictions : np.ndarray
        Signal classification probabilities, shape (n_hits,).
    targets : np.ndarray
        Reconstructed feature vectors, shape (n_hits, n_features).
    labels : np.ndarray
        True signal/noise labels (1/0), shape (n_hits,).
    threshold : float
        Classification threshold.

    Returns
    -------
    dict
        Metrics including accuracy, signal efficiency, noise rejection,
        precision, recall, F1, and reconstruction MSE.
    """
    binary_preds = (predictions > threshold).astype(np.float32)

    tp = ((binary_preds == 1) & (labels == 1)).sum()
    fp = ((binary_preds == 1) & (labels == 0)).sum()
    fn = ((binary_preds == 0) & (labels == 1)).sum()
    tn = ((binary_preds == 0) & (labels == 0)).sum()

    accuracy = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    signal_efficiency = recall
    noise_rejection = tn / (tn + fp) if (tn + fp) > 0 else 0

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "signal_efficiency": signal_efficiency,
        "noise_rejection": noise_rejection,
    }


def tracking_efficiency(predicted_tracks, true_track_ids, hit_labels,
                        min_hits=3, purity_threshold=0.5):
    """Compute tracking efficiency and fake rate.

    Parameters
    ----------
    predicted_tracks : list[list[int]]
        List of track candidates, each a list of hit indices.
    true_track_ids : np.ndarray
        True track ID for each hit (-1 for noise).
    hit_labels : np.ndarray
        Signal/noise labels (1/0) for each hit.
    min_hits : int
        Minimum hits for a valid track candidate.
    purity_threshold : float
        Minimum fraction of hits from the same true track.

    Returns
    -------
    dict
        tracking_efficiency, fake_rate, n_true_tracks, n_found, n_fake,
        duplicate_rate.
    """
    true_tracks = set(true_track_ids[true_track_ids >= 0])
    n_true = len(true_tracks)

    found_true_tracks = set()
    n_fake = 0
    n_duplicate = 0

    for candidate in predicted_tracks:
        if len(candidate) < min_hits:
            n_fake += 1
            continue

        cand_ids = true_track_ids[candidate]
        valid_ids = cand_ids[cand_ids >= 0]

        if len(valid_ids) == 0:
            n_fake += 1
            continue

        # Majority vote
        unique, counts = np.unique(valid_ids, return_counts=True)
        majority_idx = np.argmax(counts)
        majority_id = unique[majority_idx]
        purity = counts[majority_idx] / len(candidate)

        if purity >= purity_threshold:
            if majority_id in found_true_tracks:
                n_duplicate += 1
            else:
                found_true_tracks.add(majority_id)
        else:
            n_fake += 1

    n_found = len(found_true_tracks)
    eff = n_found / n_true if n_true > 0 else 0
    total_reco = len(predicted_tracks)
    fr = n_fake / total_reco if total_reco > 0 else 0
    dup = n_duplicate / total_reco if total_reco > 0 else 0

    return {
        "tracking_efficiency": eff,
        "fake_rate": fr,
        "duplicate_rate": dup,
        "n_true_tracks": n_true,
        "n_found": n_found,
        "n_fake": n_fake,
        "n_duplicate": n_duplicate,
        "n_candidates": total_reco,
    }


def fake_rate(predicted_tracks, true_track_ids, purity_threshold=0.5):
    """Compute fake track rate (convenience wrapper)."""
    labels = np.zeros(len(true_track_ids))
    labels[true_track_ids >= 0] = 1
    result = tracking_efficiency(
        predicted_tracks, true_track_ids, labels,
        purity_threshold=purity_threshold
    )
    return result["fake_rate"]
