# CaloCalibration
This repo contains scripts which take input from the Archive Tables for Source, MIP and Laser calibrations (plus others) and outputs the contents for the Reco Table.

# The Simulation

Since we do not have real data from the current system, we base our work on Simulation. A model of the calibration system is built in GEANT4 within Mu2e/Offline. Photons can be fired uniformly from this geometry to simulate the system.

To configure this simulation use the CaloCalibration/SourceCalib/fcl/Run*.fcl files.

First we run the generator, producing detector steps e.g.

```
mu2e -c CaloCalibration/SourceCalib/fcl/RunCaloCalibGun-sim.fcl --nevts=1000

```

This configures to CaloCalibGun:

```
generate: {
  module_type : CaloCalibGun
  cosmin :  -1.0
  cosmax :  1.0
  phimin :  0.0
  phimax : 2.0
  tmin  :  500.
  tmax : 100000.//Off-spill
  nDisk : 0 //Change me to run over other disk
}

```

Currently the generator must be ran one disk at a time, the last parameter chooses the disk (0 or 1).

This stage produces a file with a prefix ``dts" (detector steps). To digitize this file:

```
mu2e -c CaloCalibration/SourceCalib/fcl/RunCaloCalibGun-digi.fcl dts.owner.name.version.sequencer.art

```

where dts.owner.name.version.sequencer.art is the file from the preceeding stage. This stage produces a file prefixed ``dig" meaning it contains digitizaed calo data (as would come from the ROC during data-taking).

To analyze and perform our calibration we must run the SourceCalibAna module on this digitized ``data":

```
mu2e -c CaloCalibration/SourceCalib/fcl/RunCaloCalibGun-ana.fcl dig.owner.name.version.sequencer.art

```

This runs an analyzer which will loop over the simulated events and analyze the digis. The output is a file of N_crys histograms, these can be taken as input into the analysis. Additionally, the file can also run overlay plots of the histogrammed data of a sipm pair. This will generate an extra plot, for the start and end cry/sipm number, with normalised residuals. This is optional and will take more time to run.


# The Analysis

Now you have a file containing one histogram per crystal or sipm. The number of entries in these histograms is the number of events. To extract the calibration constants for each crystal or sim we muse run our analysis code. The analysis code is stored in "MakeAnalysisTree". This is the main script which applies the ROOFIT based maximum likelihood or chi2 fit to the crystal histograms to extract the constants.

## MakeAnalysisTree


To accumulate all the events for a given crystal the command line would look like this:
```
./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree start_crystalnumber end_crystalnumber "minimisation_method" disk number
```
To add the data comparision plots, you must add  "overlay" at the end of the command:
```
./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree start_crystalnumber end_crystalnumber "minimisation_method" disk number overlay
```
For example you can run the MakeAnalysisTree program as follows:
```
bash-5.1$  ./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 0 674 "chi2" 0
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
The arguments are the crystal ranges you want to fit to. This is followed by the method of minimization (nll or chi2) and the disk number.

## Understanding bad fits

The SourcePlotter class outputs useful histograms that can identify crystals that are significantly worse (or odd) compared to the mean and std dev of the whole calorimeter. This can be used to diagnose issues with the fits and eventually bad sipms or degraded crystals.

# The Tables

The current code outputs the results of the fit in a ROOT TTree called arXivTable.root. This will eventually be replaced with a text output in a structure that can be easily input into the archive calibration table for the sourc system.

# Combinations

In order to utilize an updated source calibration constant it must be input in a reco table and called by the digi code. This reco table uses information from both the source and cosmic tables. The exact algorithm for combining the two inputs is stored in the CaloCalibration/Combinations directory.

# Development
Current code underdevelopment by Sophie Middleton (SourceCalib+Combinations), Huma Jafree (Fitting Code) and Sam Zhou (Combinations)

