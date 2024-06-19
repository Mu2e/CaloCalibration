# CaloCalibration
This repo contains scripts which take input from the Archive Tables for Source, MIP and Laser calibrations (plus others) and outputs the contents for the Reco Table.

# The Analysis

## Module
The module SourceCalibAna should be ran on reconstructed calorimeter information. It can be ran alongside the simulation using the following:

```
mu2e -c SourceCalib/fcl/RunCaloCalibGun.fcl --nevts=2000000
```
20M events corresponds to what we expect in around 1 calibration run. It provides about 10K events per crystal.

Once the job has been run. You should run the compiled C++ script detailed below to bin by crystal and fit.

## MakeAnalysisTree

To accumulate all the events for a give crystal run the MakeAnalysisTree program as follows:

```
bash-5.1$ ./build/al9-prof-e28-p057/CaloCalibration/bin/MakeAnalysisTree chooseCrystal
========== Welcome to the Mu2e Source Calibration Analysis ==========
crystal to be analyzed (int) : 
1322
Running pre-processing .....
 Finding Hits in Crystal # 1322
 Events analyzed in this crystal : 4043
 Time take no filter crystal 24s
 Fitting Crystal # 1322
Finished pre-processing ...

```
the "chooseCrystal" option allows calibration of a single crystal e.g. 1322. If you run without that argument the code lops over all crystals, creating a .root ntuple for each.

## Fitting

Once you have the per crystal ntuples you can fit the RooFit code these....I will integrate that into the program soon.

# The Tables TODO

Our Proditions Entity, <NAME> will call upon the calibration constants in our CalEnergyCalib reco table. There is one constant per SiPM, totalling 2*1348.

The reco table values must average over all calibration modes including energy calibration from source, cosmic and laser modes.

This code currently takes the input for the source, cosmic and laser from fake datatables which are in the same format as the archive tables from real data will be.

It makes a simple linear assumption and just averages the constants. The final values are passed to the reco table.

# Development
Current code underdevelopment by Sophie Middleton 

