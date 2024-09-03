Event display to draw cosmic events on the Mu2E Calorimeter.

**To use**

On mu2e gpvms
'''
mu2einit
muse setup
/usr/bin/python3 disp_GUI.py
'''

For all other machines, use the included conda/mamba environment. Root will raise the warning:
'''
DeprecationWarning: The attribute syntax for TFile is deprecated and will be removed in ROOT 6.34. Please use TFile["sidet"] instead of TFile.sidet
'''
But this is necessary if the program needs to run on 6.30/04 on gpvms.

disp.py contains the most important calsses and a reference CLI program, while disp_GUI.py cointains the graphical user intergace elements. crystalpos.py specifies the positions of the crystals as well as the measure of their side.

At the moment, the program only considers one of the two calorimeter disks.