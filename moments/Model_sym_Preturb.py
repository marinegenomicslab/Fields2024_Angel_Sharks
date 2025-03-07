#!/bin/python3

# Example of the python file for starting a moments optimization
# This is setup for the western Gulf splitting from the Eastern Gulf before the Atlantic split from the Eastern Gulf

import sys
import os
import moments
sys.path.append('..')
import Optimize_Functions
import Angel_models
sys.path.append('.')
import Config_Make_fs_from_vcf as cfg

#Make a frequency spectrum file
print("Loading data {}".format(sys.argv[1]))
fs = moments.Spectrum.from_file(cfg.FILE)
ns = fs.sample_sizes

plabels = "West_0, West, Anc_0, Anc, East_0, East, Atl_0, Atl, m_W_Anc, m_W_E, m_E_W, T1, T2"
upper = [10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 100.0, 100.0, 100.0, 50.0, 50.0]
lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0, 0, 0, 1e-15, 1e-15]
params=None
if params is None:
    Optimize_Functions.Optimize_Routine(fs, sys.argv[1], "model_WEA", Angel_models.pop_model_sym, 1, 13, fs_folded=True, param_labels = plabels, in_upper = upper, in_lower = lower, reps = [1])
else:
    Optimize_Functions.Optimize_Routine(fs, sys.argv[1], "model_WEA", Angel_models.pop_model_sym, 1, 13, fs_folded=True, param_labels = plabels, in_upper = upper, in_lower = lower, reps = [1], in_params=params)
