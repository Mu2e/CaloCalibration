# Event display for cosmic events on the Mu2E Calorimeter.

## To use

### On mu2e gpvms
On mu2e gpvms, onnce you copied this folder and navigate to it, run:
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
## Containts
disp.py contains the most important calsses and a reference CLI program, while disp_GUI.py cointains the graphical user intergace elements. crystalpos.py specifies the positions of the crystals as well as the measure of their side.

At the moment, the program only considers one of the two calorimeter disks.
