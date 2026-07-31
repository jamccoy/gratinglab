"""
Intraclass correlation for row-group measurements

Row-group analysis produces many blaze-angle measurements per image, but they
re-measure the same physical grooves from different horizontal bands of one scan.
The rest of the statistics in this package divide by sqrt(N) over that full count,
which is only valid if the measurements are independent.

The intraclass correlation coefficient (ICC) is how that assumption gets tested:
the fraction of total variance that sits *between* row groups rather than within
them.

    ICC = 0    groups are indistinguishable; measurements behave independently
               and dividing by sqrt(N) is roughly right
    ICC = 1    all variance is between groups; the measurements within a group
               carry no new information, and the effective sample size is closer
               to the number of groups than the number of measurements

This module only measures. It changes nothing about how SEMs or p-values are
currently computed - see experimental/hierarchical_stats/ for the correction that
this measurement is meant to size.
"""
import numpy as np

# Interpretation thresholds, carried over from
# experimental/hierarchical_stats/correlation_source_analysis.py
ICC_NEGLIGIBLE = 0.1
ICC_MODERATE = 0.2


def interpret_icc(icc):
    """One-line reading of an ICC value"""
    if icc < ICC_NEGLIGIBLE:
        return "negligible - simple statistics defensible"
    if icc < ICC_MODERATE:
        return "moderate - hierarchical statistics recommended"
    return "substantial - hierarchical statistics needed"


def compute_icc(values, group_labels):
    """
    One-way ICC of `values` grouped by `group_labels`.

    Uses the same variance-ratio form as the original quick_icc_check:
    between-group variance over total variance, where within-group variance is
    the mean of the per-group variances.

    Parameters:
        values: measurements (e.g. blaze angles)
        group_labels: which group each measurement belongs to, same length

    Returns:
        dict with icc, between_var, within_var, n_groups, n_measurements,
        mean_group_size and interpretation. Returns icc = nan when there are
        fewer than two groups, since between-group variance is undefined - a
        caller must not read that as "no correlation".
    """
    values = np.asarray(values, dtype=float)
    labels = np.asarray(group_labels)

    if len(values) != len(labels):
        raise ValueError(f"values and group_labels must be the same length, "
                         f"got {len(values)} and {len(labels)}")

    unique = np.unique(labels)
    n_groups = len(unique)

    result = {
        'n_measurements': len(values),
        'n_groups': n_groups,
        'mean_group_size': len(values) / n_groups if n_groups else np.nan,
    }

    if n_groups < 2 or len(values) < 2:
        result.update({'icc': np.nan, 'between_var': np.nan,
                       'within_var': np.nan,
                       'interpretation': 'undefined - need at least 2 groups'})
        return result

    group_means = np.array([values[labels == g].mean() for g in unique])
    # Groups of size 1 contribute no within-group variance information
    group_vars = np.array([values[labels == g].var(ddof=1)
                           if (labels == g).sum() > 1 else np.nan
                           for g in unique])

    between_var = float(group_means.var(ddof=1))
    within_var = float(np.nanmean(group_vars)) if not np.all(np.isnan(group_vars)) else 0.0

    total = between_var + within_var
    icc = float(between_var / total) if total > 0 else 0.0

    result.update({
        'icc': icc,
        'between_var': between_var,
        'within_var': within_var,
        'interpretation': interpret_icc(icc),
    })
    return result


def effective_sample_size(n_measurements, mean_group_size, icc):
    """
    Sample size the measurements are actually worth, given their correlation.

    n_eff = n / (1 + (m - 1) * ICC), the standard design-effect correction, where
    m is the mean measurements per group. At ICC = 0 this returns n; as ICC rises
    it falls toward the number of groups.

    This is the number that shows how far the current SEMs are off: they are
    computed with n, and dividing by sqrt(n) instead of sqrt(n_eff) understates
    the standard error by sqrt(n / n_eff).
    """
    if not np.isfinite(icc):
        return np.nan
    design_effect = 1.0 + (mean_group_size - 1.0) * icc
    return n_measurements / design_effect if design_effect > 0 else np.nan


def sem_inflation_factor(mean_group_size, icc):
    """
    How much the reported SEM should grow: sqrt(design effect).

    A factor of 2 means the true standard error is twice what is currently
    reported, and every confidence interval is half the width it should be.
    """
    if not np.isfinite(icc):
        return np.nan
    return float(np.sqrt(1.0 + (mean_group_size - 1.0) * icc))
