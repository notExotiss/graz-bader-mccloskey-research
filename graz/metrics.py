"""Evaluation metrics for community-detection experiments."""

from __future__ import annotations

from collections import defaultdict

from sklearn.metrics import normalized_mutual_info_score


def _aligned(labels_dict, gt_dict):
    """Return (pred, truth) lists over the nodes present in BOTH dicts."""
    common = [u for u in gt_dict if u in labels_dict]
    pred = [labels_dict[u] for u in common]
    truth = [gt_dict[u] for u in common]
    return pred, truth


def nmi(labels_dict, gt_dict):
    """Normalized mutual information against ground truth."""
    pred, truth = _aligned(labels_dict, gt_dict)
    if not pred:
        return float("nan")
    return normalized_mutual_info_score(truth, pred)


def num_communities(labels_dict):
    return len(set(labels_dict.values()))


def community_f1(labels_dict, gt_dict):
    """Best-match averaged F1 between recovered and ground-truth communities.

    For every ground-truth community we find the recovered community with the
    largest F1 overlap and average those best-match F1 scores (symmetric mean of
    the two directions), which is the standard set-overlap F1 used in SNAP work.
    """
    pred, truth = _aligned(labels_dict, gt_dict)
    if not pred:
        return float("nan")
    pred_sets = defaultdict(set)
    true_sets = defaultdict(set)
    for i, (p, t) in enumerate(zip(pred, truth)):
        pred_sets[p].add(i)
        true_sets[t].add(i)

    def best_f1(a_sets, b_sets):
        total = 0.0
        for a in a_sets.values():
            best = 0.0
            for b in b_sets.values():
                inter = len(a & b)
                if inter == 0:
                    continue
                prec = inter / len(a)
                rec = inter / len(b)
                f1 = 2 * prec * rec / (prec + rec)
                if f1 > best:
                    best = f1
            total += best
        return total / len(a_sets) if a_sets else 0.0

    f_pred = best_f1(pred_sets, true_sets)
    f_true = best_f1(true_sets, pred_sets)
    return 0.5 * (f_pred + f_true)


def community_sizes(labels_dict):
    sizes = defaultdict(int)
    for c in labels_dict.values():
        sizes[c] += 1
    return sorted(sizes.values(), reverse=True)
