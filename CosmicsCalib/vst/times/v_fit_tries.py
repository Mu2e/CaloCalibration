""" Mu2E Calorimeter Calibration
    Tester for different time intervals in fitting SiPM data
    By Giacinto Boccia
    Version 1.1 | 2024-08-15
    """

import ROOT as R
import numpy as np
import tempfile
import os.path
from array import array
import concurrent.futures
from t_intervals import t_inters

#Ignoring most ROOT error messages
R.gErrorIgnoreLevel = 6001
#Load macros
R.gInterpreter.LoadMacro ('./caloreco.C')
R.gInterpreter.LoadMacro ('./AnaDriver2_1.C')

MIN_CRYSTALS = 50
results_arr = np.zeros((t_inters.shape[0], 2), dtype = np.double)
#This will store (sigma, n_couples)
result_direcory = "/exp/mu2e/data/users/gboccia/time_tests/"
#result_direcory = "./Results/"

def interval_to_name(interval : np.array) -> str:
    name = ""
    for value in interval:
        #We use "p" for positive and "m" for negative values
        if value >= 0:
            name += "p"
            name += str(abs(int(value)))
        else:
            name += "m"
            name += str(abs(int(value)))
        name += "_"
    return name[:-1]

def get_sigma(interval : np.array) -> tuple[np.double, np.double]:
    print("Trying interval", interval[0], interval[1], sep = " ")
    int_name = interval_to_name(interval)
    sigma = np.double(0)
    n_data = np.double(0)

    #Run Caloreco
    out_name_calo = result_direcory + int_name + ".root"
    if not os.path.exists(out_name_calo):
        #Skipped if output file already exists
        caloreco_obj = R.caloreco('data.list')
        caloreco_obj.Loop(out_name_calo, 0, interval[0], interval[1])
    
    #Store output name in a temp file
    tmp = tempfile.NamedTemporaryFile(dir = result_direcory, prefix = 'temp_', suffix = '.list')
    tmp.write(bytes(out_name_calo, 'utf-8'))
    tmp.flush()

    #Run AnaDriver2_1
    out_name_ana = result_direcory + int_name + "_time.root"
    if os.path.exists(out_name_ana):
        #If outupt file already exists, retrive result
        file = R.TFile.Open(out_name_ana)
        #Check if the histogram is not empty
        n_data = file["h_sig"].GetEntries()
        if n_data > 0:
            fit = file["h_sig"].GetFunction("Gauss")
            sigma = fit.GetParameter(1)
        file.Close()
    else:
        #Else run AnaDriver
        anadriver_obj = R.AnaDriver2(tmp.name)
        sigma = anadriver_obj.Loop(out_name_ana, 0, 12, 10)
        with R.TFile.Open(out_name_ana) as file:
            n_data = file.h_sig.GetEntries()
    return sigma, n_data

with concurrent.futures.ProcessPoolExecutor() as executor:
    for i, result in zip(range(t_inters.shape[0]), executor.map(get_sigma, t_inters)):
        results_arr[i] = result

#Prepare a graph
canvas = R.TCanvas("Sigma-DeltaT", "Sigma-DeltaT")
plot = R.TGraph2D()
#Fill with nonzero values that where calculated over a minimum number of crystals
n_points = 0
for i in range(t_inters.shape[0]):
    if results_arr[i, 0] > 0 and results_arr[i, 1] >= MIN_CRYSTALS:
        plot.SetPoint(n_points, t_inters[i, 0], t_inters[i, 1], results_arr[i, 0])
        n_points += 1
#Plot
plot.SetTitle("Sigma-DeltaT")
plot.GetZaxis().SetTitle("sigma")
plot.GetXaxis().SetTitle("t_start")
plot.GetYaxis().SetTitle("t_end")
plot.Draw()
canvas.SaveAs(result_direcory + "sigma_delta.root")

#Save result to a CSV file
with open (result_direcory + "time_test_reslts.csv", 'w') as file:
    print ("Start,End,Sigma,N_crystals", file = file)
    for i in range (t_inters.shape[0]):
        print (t_inters[i, 0], t_inters[i, 1], results_arr[i, 0], results_arr[i, 1], sep = ",", file = file)

input("Press any key to exit")