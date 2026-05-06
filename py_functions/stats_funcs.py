import numpy as np
from scipy.stats import ttest_ind, ttest_rel


def sigtest(yearmean1, yearmean2, timemean1, timemean2):
    """Paired t-test for fields with identical sample sizes. Masks where p > 0.1."""
    ptvals = ttest_rel(yearmean1, yearmean2, axis=0)
    diff = timemean1 - timemean2
    diff_mask = np.ma.masked_where(ptvals[1] > 0.1, diff)
    return diff, diff_mask, ptvals


def sigtest2n(yearmean1, yearmean2, timemean1, timemean2):
    """Welch's t-test for fields with unequal sample sizes. Masks where p > 0.1."""
    ptvals = ttest_ind(yearmean1, yearmean2, axis=0, equal_var=False)
    diff = timemean1 - timemean2
    diff_mask = np.ma.masked_where(ptvals[1] > 0.1, diff)
    return diff, diff_mask, ptvals


def wind_sig_mask(diff_mask_u, diff_mask_v):
    """Return (u_sig, v_sig) arrays masked to grid points where BOTH components are significant.

    Parameters
    ----------
    diff_mask_u, diff_mask_v : np.ma.MaskedArray
        Significance-masked u and v wind difference fields from sigtest / sigtest2n.
    """
    both_sig = (diff_mask_u.mask == False) & (diff_mask_v.mask == False)
    u_sig = np.where(both_sig, diff_mask_u.data, np.nan)
    v_sig = np.where(both_sig, diff_mask_v.data, np.nan)
    return u_sig, v_sig
