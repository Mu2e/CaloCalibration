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

The code assumes two input archive tables:

* source ( stored as Source.txt)  - the source calibration
* cosmic ( stored as Cosmics.txt) - the cosmic/mip calibration

Originally it was assumed these would have the same structure. We need to edit the code for the source.

## Producing a reco table

Within ```CaloCalibTableMaker.hh``` we have two structs defined: the ArchiveTable and RecoTable.

The point is that the ArchiveTable is channels specific and the RecoTable is the simple combined output.

The main code can be configured (assuming Cosmic.txt and Source.txt are local):

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker 1 1
```

For just cosmics:

```
./build/al9-prof-e29-p082/CaloCalibration/bin/CaloCalibTableMaker 1
```

## Development

Currently being developed by Sam Zhou and Sophie Middleton
