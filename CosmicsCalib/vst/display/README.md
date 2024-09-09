# Cosmic Event Display for the Mu2E Calorimeter.

## To run

### On mu2e gpvms
On mu2e gpvms, once you copied this folder and navigate to it, run:
```
mu2einit
muse setup
/usr/bin/python3 disp_GUI.py
```
### Other machines
For all other machines, use the included **conda/mamba environment**. Root will then raise the warning:
```
DeprecationWarning: The attribute syntax for TFile is deprecated and will be removed in ROOT 6.34.
Please use TFile["sidet"] instead of TFile.sidet
```
It should not break the program, the use of this command is necessary if the program needs to run on Root 6.30/04 on gpvms.

## How to use
The programs propts the use to select a file to display, this file must contain a Root Tree named "sidet" with reconstructed events. Then the user sets cuts on the reconstructed Q to be considered in fits, the number of hits, the Chi squared of a linear fit and the policy about vertical tracks that it wants to adopt. When dealing wih vertical tracks the Chi squared cut is ignored.

Each time that the Go button is pressed, the program will try to find an event that meets the cryteria with a number greather than the current (displayd in the run and event fields), if tree ends and before a match, the program will end the search without sowing a graph, but you can still use it by typing 0 in the run and event fields and pressing Go twice. Each time that ueser types a value in the run and vent fields, the subsequent use of Go will ignore all the exhisting cuts and try to display that user selected event (or the closest if it doesn't exhist).

## Contains
disp.py contains the most important classes and a reference CLI program, while disp_GUI.py contains the graphical user interface elements. crystalpos.py specifies the positions of the crystals as well as the measure of their side.

At the moment, the program only considers one of the two calorimeter disks.
