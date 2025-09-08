# CaloCalibration
This repo contains scripts which take input from the Archive Tables for Source, MIP and Laser calibrations (plus others) and outputs the contents for the Reco Table.

# The Simulation

Since we do not have real data from the current system, we base our work on Simulation. A model of the calibration system is built in GEANT4 within Mu2e/Offline. Photons can be fired uniformly from this geometry to simulate the system.

To configure this simulation use the CaloCalibration/SourceCalib/fcl/Run*.fcl files.

```
mu2e -c CaloCalibration/SourceCalib/fcl/RunCaloCalibGunDigi.fcl
```

This runs an analyzer which will loop over the simulated events and analyze the digis. The output is a file of N_crys histograms, these can be taken as input into the analysis. Additionally, the file can also run overlay plots of the histogrammed data of a sipm pair. This will generate an extra plot, for the start and end cry/sipm number, with normalised residuals. This is optional and will take more time to run.

# The Analysis

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
The arguments are the crystal ranges you want to fit to.

## Fitting

Once you have the per crystal ntuples you can fit the RooFit code these. The SourceFitter class utilizes RooFit to fit an individual crystal. The user can chose between an nll or chi2 fit by setting the  third input arg to either "nll" or "chi2".

The code is currently setup to loop over the chosen crystals, fitting to each one and storing the resulting fit parameters. An intermediate .root file is created. In theory, this would be our archive table.

We then use the SourcePlotter class to plot the overall distribution of constants and other features over all the crystals fitted. This allows us to look for trends and possible sources of bad fits.

# The Tables TODO

Our Proditions Entity, <NAME> will call upon the calibration constants in our CalEnergyCalib reco table. There is one constant per SiPM, totalling 2*1348.

The reco table values must average over all calibration modes including energy calibration from source, cosmic and laser modes.

This code currently takes the input for the source, cosmic and laser from fake datatables which are in the same format as the archive tables from real data will be.

It makes a simple linear assumption and just averages the constants. The final values are passed to the reco table.

# Development
Current code underdevelopment by Sophie Middleton and Huma Jafree

