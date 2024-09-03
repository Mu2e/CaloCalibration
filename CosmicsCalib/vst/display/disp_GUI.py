""" Mu2E Calorimeter Calibration
    Track event display GUI version 1.1
    by Giacinto Boccia
    2024-00-01
    """

import tkinter as tk
from tkinter import filedialog
import ROOT as R
import disp

def control_panel() -> None:
    window = tk.Tk()
    window.title("Event Display")
    par_fields_arr = [] 
    parameters = dict([ ('Run Number', 0),
                        ('Event Number', 0), 
                        ('Q Threshold', 4000),
                        ('Minimum Hits', 6), 
                        ('Maximum ChiSq', 10)])
    vert_options : tuple[str] = ("Include vertical tracks", "Exclude vertical tracks", "Only vertical tracks")
    selected_vertical = tk.StringVar(window)
    selected_vertical.set('Include vertical tracks')
    
    #Ask for the file to open
    file_path = filedialog.askopenfilename()
    file = R.TFile.Open(file_path)
    tree = file["sidet"]
    tree.BuildIndex("nrun","evnum")
    calo = disp.Disk(0)

    for i, (name, value) in enumerate(parameters.items()):
        tk.Label(window, text=name).grid(row=i, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i, column=1)
        entry.insert(0, value)
        par_fields_arr.append(entry)
    v_mode_menu = tk.OptionMenu(window, selected_vertical, *vert_options).grid(row=4, column=0)

    def go_action() -> None:
        #Collect parameters
        for field, p_name in zip(par_fields_arr, parameters):
            if p_name == 'Event Number' or p_name == 'Run Number' or p_name == 'Minimum Hits':
                parameters[p_name] = int(field.get())
            else:
                parameters[p_name] = float(field.get())
        match selected_vertical.get():
            case "Include vertical tracks":
                vert_mode = 'i'
            case "Exclude vertical tracks":
                vert_mode = 'e'
            case "Only vertical tracks":
                vert_mode = 'o'
            
        #Find next event to display
        entry_n = tree.GetEntryNumberWithBestIndex(parameters["Run Number"], parameters["Event Number"])
        while entry_n < tree.GetEntries() - 1:
            calo.empty()
            entry_n += 1
            if disp.single_event_q(calo, tree, 
                                   entry_n, 
                                   parameters['Q Threshold'], 
                                   parameters['Minimum Hits'], 
                                   parameters['Maximum ChiSq'],
                                   vert_mode):
                break
        tree.GetEntry(entry_n)
        parameters["Run Number"] = tree.nrun
        parameters["Event Number"] = tree.evnum
        par_fields_arr[0].delete(0, tk.END)
        par_fields_arr[0].insert(0, parameters["Run Number"])
        par_fields_arr[1].delete(0, tk.END)
        par_fields_arr[1].insert(0, parameters["Event Number"])

    def terminate_action() -> None:
        print("Termination button pressed")
        window.destroy()
        
    def average_action() -> None:
        averages_control_panel(tree)
        
    def td_action() -> None:
        disp.signle_event_td(calo, tree, parameters["Event Number"])

    tk.Button(window, text="Go",        command=go_action)          .grid(row=6, column=0)
    tk.Button(window, text="T Diff.",   command=td_action)          .grid(row=6, column=1)
    tk.Button(window, text="Terminate", command=terminate_action)   .grid(row=7, column=0)
    tk.Button(window, text="Averages",  command=average_action)     .grid(row=7, column=1)

    window.mainloop()
    
def averages_control_panel(tree) -> None:
    window = tk.Tk()
    window.title("Averages Display")
    calo = disp.Disk(0)
    disp.load_tree(calo, tree)
    
    def terminate_action() -> None:
        window.destroy()
        
    def td_action() -> None:
        calo.draw_tdif(plot_name= "Average Time Differences")
        
    def mean_q_action() ->None:
        calo.draw_q(fits= False, plot_name= "Agerage Q Values")
        
    def num_action() -> None:
        calo.draw_hitcount()
        
    
    tk.Button(window, text="Time Differences",  command=td_action)          .grid(row=0, column=0)
    tk.Button(window, text="Number of Hits",    command=num_action)         .grid(row=0, column=1)
    tk.Button(window, text="Close Averages",    command=terminate_action)   .grid(row=1, column=0)
    tk.Button(window, text="Average Qs",        command=mean_q_action)      .grid(row=1, column=1)

if __name__ == "__main__":
    control_panel()
