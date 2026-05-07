import numpy as np
from scipy.stats import ttest_ind, ttest_rel


def sigtest(yearmean1, yearmean2, timemean1, timemean2):
    """Paired t-test for fields with identical sample sizes. Returns NaN where p > 0.1."""
    ptvals = ttest_rel(yearmean1, yearmean2, axis=0)
    diff = timemean1 - timemean2
    diff_mask = diff.where(ptvals[1] <= 0.1)
    return diff, diff_mask, ptvals


def sigtest2n(yearmean1, yearmean2, timemean1, timemean2):
    """Welch's t-test for fields with unequal sample sizes. Returns NaN where p > 0.1."""
    ptvals = ttest_ind(yearmean1, yearmean2, axis=0, equal_var=False)
    diff = timemean1 - timemean2
    diff_mask = diff.where(ptvals[1] <= 0.1)
    return diff, diff_mask, ptvals


def wind_sig_mask(diff_mask_u, diff_mask_v):
    """Return (u_sig, v_sig) DataArrays masked to grid points where BOTH components are significant.

    Parameters
    ----------
    diff_mask_u, diff_mask_v : xr.DataArray
        Significance-masked u and v wind difference fields from sigtest / sigtest2n.
        NaN where not significant.
    """
    both_sig = diff_mask_u.notnull() & diff_mask_v.notnull()
    return diff_mask_u.where(both_sig), diff_mask_v.where(both_sig)
