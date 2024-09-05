""" Time Calibration of the Mu2e Calorimeter
    by Giacinto boccia
    version 0.1 | 2024-09-05
"""
import numpy as np
import dask as dk
import uproot
import quantities as pq

hits_path = input("Hits file to process:")
cal_path = input("Starting caibration file:") or False

