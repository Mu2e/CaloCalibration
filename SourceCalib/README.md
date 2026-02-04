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

Now you have a file containing one histogram per crystal or sipm. The number of entries in these histograms is the number of events. To extract the calibration constants for each crystal or sipm we must run our analysis code. The analysis code is stored in "MakeAnalysisTree". This is the main script which applies the ROOFIT based maximum likelihood or chi2 fit to the crystal histograms to extract the constants.

## MakeAnalysisTree
The MakeAnalysisTree executable performs source calibration analysis. It can either fit crystal energy spectra (Data mode) or extract truth information (MC mode).
### Usage
The command requires 4 mandatory arguments followed by any optional flags.
### Arguments

| Argument    | Type   | Description                                                                 |
| :---        | :---   | :---                                                                        |
| **start_cry** | int    | The starting crystal ID to analyze.                                         |
| **end_cry** | int    | The final crystal ID (exclusive limit).                                     |
| **alg** | string | The fitting algorithm: "nll" (Log-Likelihood) or "chi2" (Chi-Square).       |
| **disk** | int    | The disk number (0 or 1).                                                   |
*Note: Currently the disk can be toggled between sipms (disk 0 or 1) vs crystals (disk 2 or 3-- representing 0 and 1 respectively) as input. This is a temporary switch being used for studies and will be removed before the final data run.*
### Optional Flags
Add one or more of these flags at the end of the command (order does not matter):

* `overlay`: Generates plots comparing the Data histogram vs. the Fit function.
* `contour`: Generates 2D parameter scan plots (e.g., Peak vs. Width) to check correlations.
* `mc`: Runs **MC Truth Mode**.
    * *Note: When `mc` is used, the program skips the fitting process entirely. You must still provide a dummy value for the `<alg>` argument (e.g., "chi2") to satisfy the command structure.*
    
To accumulate all the events for a given crystal the command line would look like this:
```
./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree <start_cry> <end_cry> <alg> <disk> [flags]
```
### Examples

**1. Standard Analysis (Chi2 Fit)**
Accumulate and fit events for crystals 0 to 674 on disk 0 using Chi-Square minimization:

./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 0 674 "chi2" 0

**2. Analysis with Overlay Plots**
Perform an NLL fit and save the data comparison plots:

./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 1322 1323 "nll" 0 overlay
**3. Analysis with Contour Plots**
Perform a fit and generate 2D error contours (useful for debugging parameter correlations):

./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 100 101 "nll" 0 contour

**4. MC Truth Extraction**
Run in MC mode to generate the truthinfo.root file. The fitting algorithm is ignored, but must be present:

./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 0 674 "chi2" 0 mc

**5. Combined Flags**
Run NLL fit, generate overlays, AND generate contours:

./build/al9-prof-e28-p056/CaloCalibration/bin/MakeAnalysisTree 100 105 "nll" 0 overlay contour

*Example Output:*
========== Welcome to the Mu2e Source Calibration Analysis ==========
crystal to be analyzed (int) : 1322
Running pre-processing .....
 Finding Hits in Crystal # 1322
 Events analyzed in this crystal : 4043
 Fitting Crystal # 1322
Finished pre-processing ...


## Understanding bad fits

The SourcePlotter class outputs useful histograms that can identify crystals that are significantly worse (or odd) compared to the mean and std dev of the whole calorimeter. This can be used to diagnose issues with the fits and eventually bad sipms or degraded crystals.

# The Tables

The current code outputs the results of the fit in a ROOT TTree called arXivTable.root. This will eventually be replaced with a text output in a structure that can be easily input into the archive calibration table for the sourc system.

# Combinations

In order to utilize an updated source calibration constant it must be input in a reco table and called by the digi code. This reco table uses information from both the source and cosmic tables. The exact algorithm for combining the two inputs is stored in the CaloCalibration/Combinations directory.

# Development
Current code underdevelopment by Sophie Middleton (SourceCalib+Combinations), Huma Jafree (Fitting Code) and Sam Zhou (Combinations)

