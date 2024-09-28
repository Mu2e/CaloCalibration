# Cosmic Event Calibration for the Mu2E Calorimeter.

## To run

### In general
Use the included **conda/mamba environment** it should come with all the dependencies, and use t_cal.py, if you can avoid the graphs option the program gets faster. 

### On the EAF
Sadly, unsolved issues prevent the parelellization of some of the operations on the EAF, therefore t_cal.ipynb is provided. It doesn't use ProcessPoolExecutor on those operations, therefore it runs on the EAF, the untallelelized operations are not in the main loop but they are called once on all the events (one is mandatory, the second only occurs if graphs are selected).
