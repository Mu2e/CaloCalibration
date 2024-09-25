""" Verifier for t_val in Mu2e calorimeter VST
    by Giacinto Boccia
    version 0.1 | 2024-09-25 """
    
import matplotlib.pyplot as plt
import uproot
import awkward as ak
from concurrent.futures import ProcessPoolExecutor

def get_t_val_interval(event) -> float | None:
    t_vals = []
    for channel in chs:
        for hit_r, hit_c, hit_s, hit_t in zip(event.iRow, event.iCol, event.SiPM, event.Tval):
            if (hit_r, hit_c, hit_s) == channel:
                t_vals.append(hit_t)
                
    if len(t_vals) == 2:
        return abs(t_vals[0] - t_vals[1])

path = input("File (.root) to open: ") + ":sidet"
chs = []
chs.append(tuple(int(ind) for ind in input("First channel (row col sipm): ").split()))
chs.append(tuple(int(ind) for ind in input("Second channel (row col sipm): ").split()))

branches = ("iRow", "iCol", "SiPM", "Tval")
with uproot.open(path) as file:
    tree = file.arrays(filter_name = branches)
      
with ProcessPoolExecutor() as executor:
    intervals = tuple(executor.map(get_t_val_interval, tree))
intervals = tuple(filter(lambda interval: interval is not None, intervals))
plt.ion()
plt.hist(intervals, bins = 1000, range= (0, 1000))
plt.title("Tval differneces between " + str(chs[0]) + " and " + str(chs[1]))

input("Press any key to exit")