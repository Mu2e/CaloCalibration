# Combinations

This directory contains the code to combine calibration from the Cosmic and Source origins into one calibration constant.

## Concept

Both Cosmics and Source calibrations will take place independently. Source has access to 6.13MeV photons, Cosmics uses 20MeV cosmics.

We can compare the two constants to understand how things change with energy.

## `CaloCalibTableMaker`

The current code is written in C++, it is compiled and built when the user runs:

```
muse build -j N
```

where N = number of cores.

If the build completes successfully than the exectuable will appear in the ```build``` repo. It can be run as follows:

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker
```
where ```al9-prof-e29-p082``` is an example of an environment build.

To acutally use the code, you will need the latest tables locally. Use openTool:

```
openTool get-calibration --name CalCosmicEnergyCalib --run 106521 &> Cosmics.txt
```

this is an example of how to extract the most up to date cosmics table.

The code now applies the calibration-combination algorithm:

1. Fit quality gate (`p > 0.01`) for cosmic and source fits.
2. Case split:
   - neither valid: fallback to nominal (`R0`)
   - one valid: use that estimator
   - both valid: inverse-variance weighted average
3. Compatibility test for the two-method case (`p_compat > 0.05`), otherwise fallback to nominal.

Expected inputs:

* `Cosmics.txt` (`CalCosmicEnergyCalib`)
* `Source.txt` (`CalSourceEnergyCalib`)
* `nominal.txt` (`CalEnergyCalib`, used for fallback `R0`)

## Producing a reco table

Within `CaloCalibTableMaker.hh` the key row structs are:

* `ArchiveFitRow` for per-channel fit inputs (cosmic/source peak, error, chi2, ndf)
* `CalEnergyCalibRow` for nominal fallback constants (`R0`)
* `CombinedCalibRow` for the final combined output with uncertainty and status

Default run (both methods enabled, default filenames above):

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker 1 1
```

Run with explicit paths:

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker 1 1 Cosmics.txt Source.txt nominal.txt RecoTable.txt CalEnergyCalibInfo.txt
```

For just cosmics (source disabled):

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker 1 0
```

Outputs:

* `RecoTable.txt` with `TABLE CalEnergyCalib`
* `CalEnergyCalibInfo.txt` with status/error metadata per channel

## Development

Currently being developed by Sam Zhou and Sophie Middleton
