""" Mu2E Calorimeter Calibration
    Track event display GUI version 1.0
    by Giacinto Boccia
    2024-08-31
    """

import tkinter as tk
from tkinter import filedialog
import ROOT as R
import disp

def control_panel() -> None:
    window = tk.Tk()
    window.title("Event Display")
    par_fields_arr = []
    parameters : dict = dict([('Event number', 0), 
                              ('Q threshold', 4000), 
                              ('Minimum hits', 6), 
                              ('Maximum ChiSq', 10)])
    
    
    #Ask for the file to open
    file_path = filedialog.askopenfilename()
    file = R.TFile.Open(file_name)
    tree = file.sidet
    calo = disp.calo

    for i, (name, value) in enumerate(parameters.items()):
        tk.Label(window, text=name).grid(row=i, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i, column=1)
        entry.insert(0, value)
        par_fields_arr.append(entry)

    def go_action() -> None:
        #Collect parameters
        for field, p_name in zip(par_fields_arr, parameters):
            parameters[p_name] = field.get()
            
        while parameters['Event Number'] < tree.GetEntries():
            calo.empty()
            if disp.single_event_q(calo, tree, parameters["Event Number"], parameters['Q threshold'], parameters['Minimum hits'], parameters['Maximum ChiSq']):
                break
            else:
                parameters["Event Number"] += 1        

    def terminate_action() -> None:
        print("Termination button pressed")
        window.destroy()
        
    def average_action() -> None:
        averages_control_panel(file_path)
        window.destroy()
        
    def td_action() -> None:
        pass

    tk.Button(window, text="Go", command=go_action).grid(row=4, column=0)
    tk.Button(window, text="T Diff.", command=td_action).grid(row=4, column=1)
    tk.Button(window, text="Terminate", command=terminate_action).grid(row=5, column=0)
    tk.Button(window, text="Averages", command=average_action).grid(row=5, column=1)

    window.mainloop()
    
def averages_control_panel(file_path) -> None:
    window = tk.Tk()
    window.title("Averages Display")
    
    def terminate_action() -> None:
        print("Termination button pressed")
        window.destroy()
        
    def td_action() -> None:
        print("td button pressed")
        
        
    def num_action() -> None:
        print("num button pressed")
        
    
    tk.Button(window, text = "Time Differences", command = td_action).grid(row = 0, column = 0)
    tk.Button(window, text = "Number of Hits", command = num_action).grid(row = 0, column = 1)
    tk.Button(window, text = "Terminate", command = terminate_action).grid(row = 5, column = 0)

if __name__ == "__main__":
    control_panel()
