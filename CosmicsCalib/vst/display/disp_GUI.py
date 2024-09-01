""" Mu2E Calorimeter Calibration
    Track event display GUI version 1.0
    by Giacinto Boccia
    2024-00-01
    """

import tkinter as tk
from tkinter import filedialog
import ROOT as R
import disp

parameters = dict([ ('Event Number', 0), 
                    ('Q threshold', 4000), 
                    ('Minimum hits', 6), 
                    ('Maximum ChiSq', 10)])

def control_panel() -> None:
    window = tk.Tk()
    window.title("Event Display")
    par_fields_arr = []   
    
    #Ask for the file to open
    file_path = filedialog.askopenfilename()
    file = R.TFile.Open(file_path)
    tree = file["sidet"]
    calo = disp.Disk(0)

    for i, (name, value) in enumerate(parameters.items()):
        tk.Label(window, text=name).grid(row=i, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i, column=1)
        entry.insert(0, value)
        par_fields_arr.append(entry)

    def go_action() -> None:
        #Collect parameters
        for field, p_name in zip(par_fields_arr, parameters):
            if p_name == 'Event Number' or p_name == 'Minimum hits':
                parameters[p_name] = int(field.get())
            else:
                parameters[p_name] = float(field.get())
            
        while parameters['Event Number'] < tree.GetEntries() - 1:
            calo.empty()
            parameters["Event Number"] += 1
            if disp.single_event_q(calo, tree, 
                                   parameters["Event Number"], 
                                   parameters['Q threshold'], 
                                   parameters['Minimum hits'], 
                                   parameters['Maximum ChiSq']):
                par_fields_arr[0].delete(0, tk.END)
                par_fields_arr[0].insert(0, parameters["Event Number"])
                break      

    def terminate_action() -> None:
        print("Termination button pressed")
        window.destroy()
        
    def average_action() -> None:
        averages_control_panel(tree)
        
    def td_action() -> None:
        disp.signle_event_td(calo, tree, parameters["Event Number"])

    tk.Button(window, text="Go",        command=go_action)          .grid(row=4, column=0)
    tk.Button(window, text="T Diff.",   command=td_action)          .grid(row=4, column=1)
    tk.Button(window, text="Terminate", command=terminate_action)   .grid(row=5, column=0)
    tk.Button(window, text="Averages",  command=average_action)     .grid(row=5, column=1)

    window.mainloop()
    
def averages_control_panel(tree) -> None:
    window = tk.Tk()
    window.title("Averages Display")
    calo = disp.Disk(0)
    disp.load_tree(calo, tree)
    
    def terminate_action() -> None:
        window.destroy()
        
    def td_action() -> None:
        calo.draw_tdif("Average Time Differences")
        
        
    def num_action() -> None:
        calo.draw_hitcount()
        
    
    tk.Button(window, text="Time Differences",  command=td_action)          .grid(row=0, column=0)
    tk.Button(window, text="Number of Hits",    command=num_action)         .grid(row=0, column=1)
    tk.Button(window, text="Close Averages",    command=terminate_action)   .grid(row=5, column=0)


if __name__ == "__main__":
    control_panel()
