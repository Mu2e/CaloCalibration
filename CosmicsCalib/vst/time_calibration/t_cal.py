""" Time Calibration of the Mu2e Calorimeter
    by Giacinto boccia
    version 0.1 | 2024-09-05
"""
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import uproot
import quantities as pq
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

#Cut parameters
HITNUM_CUT = 6
V_MIN_CUT = 338
V_MAX_CUT = 675
COS_THETA_CUT = 0.2
CHI_ON_NDF_CUT = 2
#constants
CORR_FACTOR = 0.5
N_ROWS = 36
N_COLUMNS = 28
N_SIPMS = 2

#Input
hits_path = input("Hits file to process:")
cal_path = input("Starting caibration file:") or False
n_runs = int(input("Iterations to perform:"))
save_f_name = input("Calibration file to save [Deafault <hits>_t_calibration.csv]:") or False
if not save_f_name:
    save_f_name = hits_path[ : -5] + "_t_calibration.csv"
#Name of the tree inside the file
hits_path += ":sidet"