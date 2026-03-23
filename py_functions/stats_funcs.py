import os
import sys
import xarray as xr
import netCDF4 as nc
import numpy as np
import metpy.calc as mp
from scipy.stats import ttest_ind, ttest_rel
from datetime import datetime

#############################

def sigtest(yearmean1,yearmean2,timemean1,timemean2):
	"""
	Determine statistical significance of two fields with identical sample sizes
	"""
	ptvals = ttest_rel(yearmean1,yearmean2, axis=0)
	diff = timemean1-timemean2
	diff_mask = np.ma.masked_where(ptvals[1] > 0.1,diff)
	return diff, diff_mask, ptvals

def sigtest2n(yearmean1,yearmean2,timemean1,timemean2):
	"""
	Determine statistical significance of two fields with unequal sample sizes
	"""
	ptvals = ttest_ind(yearmean1,yearmean2, axis=0, equal_var = False)
	diff = timemean1-timemean2
	diff_mask = np.ma.masked_where(ptvals[1] > 0.1,diff)
	return diff, diff_mask, ptvals

