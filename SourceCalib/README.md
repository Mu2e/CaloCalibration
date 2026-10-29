# CaloCalibration
This repo contains scripts which take input from the Archive Tables for Source, MIP and Laser calibrations (plus others) and outputs the contents for the Reco Table.

# The Simulation

Since we do not have real data from the current system, we base our work on Simulation. A model of the calibration system is built in GEANT4 within Mu2e/Offline. Photons can be fired uniformly from this geometry to simulate the system.

To configure this simulation use the CaloCalibration/SourceCalib/simulation/Run*.fcl files.

First we run the generator, producing detector steps e.g.

```
mu2e -c CaloCalibration/SourceCalib/simulation/RunCaloCalibGun-sim.fcl --nevts=1000

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
mu2e -c CaloCalibration/SourceCalib/simulation/RunCaloCalibGun-dig.fcl dts.owner.name.version.sequencer.art

```

where dts.owner.name.version.sequencer.art is the file from the preceeding stage. This stage produces a file prefixed ``dig" meaning it contains digitizaed calo data (as would come from the ROC during data-taking).

To analyze and perform our calibration we must run the SourceCalibAna module on this digitized ``data":

```
mu2e -c CaloCalibration/SourceCalib/simulation/RunCaloCalibGun-ana.fcl dig.owner.name.version.sequencer.art

```

This runs an analyzer which will loop over the simulated events and analyze the digis. The output is a file of N_crys histograms, these can be taken as input into the analysis. Additionally, the file can also run overlay plots of the histogrammed data of a sipm pair. This will generate an extra plot, for the start and end cry/sipm number, with normalised residuals. This is optional and will take more time to run.


# The Analysis

Now you have a file containing one histogram per crystal or sipm. The number of entries in these histograms is the number of events. To extract the calibration constants for each crystal or sipm we must run our analysis code. The analysis code is stored in "MakeAnalysisTree". This is the main script which applies the ROOFIT based chi2 fit to the crystal histograms to extract the constants.

## MakeAnalysisTree
The MakeAnalysisTree executable performs the source calibration analysis. It either fits the ADC spectrum of each SiPM or crystal (fit mode), or extracts truth information from the MC truth histograms (MC mode).

The input file is set by `filepath` at the top of `src/MakeAnalysisTree_main.cc`. To run on a different file, change it there and rebuild.

### Usage
```
./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree <start_id> <end_id> <disk> <readout> [flags]
```
All four arguments are required. The fit is always a χ² fit: the earlier `<alg>` argument and its `nll` (log-likelihood) option were removed because each nll fit took about 100× longer than chi2.

### Arguments
| Argument | Type | Description |
| :--- | :--- | :--- |
| **start_id** | int | First SiPM or crystal ID to analyse. |
| **end_id** | int | Last ID (exclusive). |
| **disk** | int | Disk number (0 or 1). Required, but it doesn't currently affect the analysis: the input file contains both disks, so the disk is chosen by the ID range (crystals 0–673 / SiPMs 0–1347 are disk 0). |
| **readout** | string | `sipm`: fit each SiPM histogram (IDs 0–2695). `crystal`: fit each crystal histogram, the sum of its SiPM pair (IDs 0–1347). |

### Optional flags
Add any of these after the four arguments, in any order:

* `singlesided`: fit with a single-sided Crystal Ball (default is double-sided).
* `overlay`: overlay each even/odd SiPM pair in the range, with a normalised residual panel (`mu2e_simu_hist_<i>_<i+1>.root`). Use with `sipm` readout and a range of at least two IDs.
* `contour <x> <y>`: 2D χ² scan of two fit parameters for each fitted ID. Allowed: `peak`, `n_full`, `n_1st`, `n_2nd`, `n_bkg`.
* `mc`: MC truth mode. Reads the MC truth energy histograms (`crystals_edep_truth`, in MeV) instead of the ADC histograms, and skips the fit. IDs are crystal IDs (0–1347) and `readout` is ignored, but it must still be given. Writes `truthinfo.root`.

If `truthinfo.root` is in the working directory during a fit run, each fit is compared with the MC truth for its crystal and the `truth_*` branches of `arXivTable.root` are filled.

### Outputs (fit mode)
* `arXivTable.root`: TTree `covar`, one entry per fitted ID with all fit parameters.
* `mu2e_simu_fitSpec_chi2_<CsI|LYSO>_<id>.root`: fit plot for each ID.
* `globalPlots.root` and the SourcePlotter summary plots.
* `Fit_Summary.txt`: how many IDs converged at each fit stage, bad χ², and convergence failures.

### Examples
**1. Chi2 fit of the disk-0 crystals**

    ./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree 0 674 0 crystal

**2. SiPM fit with even/odd overlay plots**

    ./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree 1322 1324 0 sipm overlay

**3. Fit with a contour scan of peak vs. full-peak yield**

    ./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree 100 101 0 sipm contour peak n_full

**4. MC truth extraction (writes truthinfo.root)**

    ./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree 0 674 0 crystal mc

**5. Single-sided fit with overlays and a contour scan**

    ./build/<build-dir>/CaloCalibration/bin/MakeAnalysisTree 100 104 0 sipm singlesided overlay contour peak n_bkg

## Understanding bad fits

The SourcePlotter class outputs useful histograms that can identify crystals that are significantly worse (or odd) compared to the mean and std dev of the whole calorimeter. This can be used to diagnose issues with the fits and eventually bad sipms or degraded crystals.

# The Tables

The current code outputs the results of the fit in a ROOT TTree called arXivTable.root. This will eventually be replaced with a text output in a structure that can be easily input into the archive calibration table for the sourc system.

# Combinations

In order to utilize an updated source calibration constant it must be input in a reco table and called by the digi code. This reco table uses information from both the source and cosmic tables. The exact algorithm for combining the two inputs is stored in the CaloCalibration/Combinations directory.

# Development
Current code underdevelopment by Sophie Middleton (SourceCalib+Combinations), Huma Jafree (Fitting Code) and Sam Zhou (Combinations)

