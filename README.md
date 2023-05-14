# CaloCalibration
This repo contains scripts which take input from the Archive Tables for Source, MIP and Laser calibrations (plus others) and outputs the contents for the Reco Table.

# The Idea

Our Proditions Entity, <NAME> will call upon the calibration constants in our CalEnergyCalib reco table. There is one constant per SiPM, totalling 2*1348.

The reco table values must average over all calibration modes including energy calibration from source, cosmic and laser modes.

This code currently takes the input for the source, cosmic and laser from fake datatables which are in the same format as the archive tables from real data will be.

It makes a simple linear assumption and just averages the constants. The final values are passed to the reco table.

# Development
Current code underdevelopment by Sophie Middleton 

