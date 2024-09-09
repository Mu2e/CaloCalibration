""" Mu2E Calorimeter Calibration
    Track event display GUI version 1.2
    by Giacinto Boccia
    2024-00-03
    """

import tkinter as tk
from tkinter import filedialog
import ROOT as R
import disp

def control_panel() -> None:
    window = tk.Tk()
    window.title("Event Display")
    par_fields_arr = [] 
    parameters = dict([('Run Number',      0),
                       ('Event Number',    0), 
                       ('Q Threshold',     4000.),
                       ('Minimum Hits',    6), 
                       ('Maximum ChiSq',   10.)])
    vert_options = dict([("Include vertical tracks", "i"),
                         ("Exclude vertical tracks", "e"),
                         ("Only vertical tracks",    "o")])
    selected_vertical = tk.StringVar(window)
    selected_vertical.set('Include vertical tracks')
    
    #Open file
    file_path = filedialog.askopenfilename()
    file = R.TFile.Open(file_path)
    tree = file.sidet
    tree.BuildIndex("nrun","evnum")

    #Declare fields
    for i, (name, value) in enumerate(parameters.items()):
        tk.Label(window, text=name).grid(row=i, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i, column=1)
        entry.insert(0, str(value))
        par_fields_arr.append(entry)
    v_mode_menu = tk.OptionMenu(window, selected_vertical, *vert_options).grid(row=5, column=0)

    def go_action() -> None:
        #Collect parameters
        user_set_ev : bool = False
        for field, p_name in zip(par_fields_arr, parameters):
            if p_name == 'Event Number' or p_name == 'Run Number' or p_name == 'Minimum Hits':
                if p_name == 'Event Number' or p_name == 'Run Number':
                    if parameters[p_name] != int(field.get()):
                        user_set_ev = True
                parameters[p_name] = int(field.get())
            else:
                parameters[p_name] = float(field.get())
        vert_mode = vert_options[selected_vertical.get()]
            
        #Display the next event
        parameters["Run Number"], parameters ["Event Number"] = disp.single_event_q(tree, 
                                                                                    parameters["Run Number"],
                                                                                    parameters["Event Number"],
                                                                                    parameters["Q Threshold"],
                                                                                    parameters["Minimum Hits"],
                                                                                    parameters["Maximum ChiSq"],
                                                                                    vert_mode,
                                                                                    user_set_ev)
        
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
        disp.signle_event_td(tree, parameters["Run Number"], parameters["Event Number"])

    tk.Button(window, text="Go",        command=go_action)          .grid(row=6, column=0)
    tk.Button(window, text="T Diff.",   command=td_action)          .grid(row=6, column=1)
    tk.Button(window, text="Terminate", command=terminate_action)   .grid(row=7, column=0)
    tk.Button(window, text="Averages",  command=average_action)     .grid(row=7, column=1)

    window.mainloop()
    
def averages_control_panel(tree : R.TTree) -> None:
    window = tk.Tk()
    window.title("Averages Display")
    disp.load_tree(tree)
    
    def terminate_action() -> None:
        window.destroy()
        
    def td_action() -> None:
        disp.calo.draw_tdif(plot_name= "Average Time Differences")
        
    def mean_q_action() ->None:
        disp.calo.draw_q(fits= False, plot_name= "Agerage Q Values")
        
    def num_action() -> None:
        disp.calo.draw_hitcount()
        
    
    tk.Button(window, text="Time Differences",  command=td_action)          .grid(row=0, column=0)
    tk.Button(window, text="Number of Hits",    command=num_action)         .grid(row=0, column=1)
    tk.Button(window, text="Close Averages",    command=terminate_action)   .grid(row=1, column=0)
    tk.Button(window, text="Average Qs",        command=mean_q_action)      .grid(row=1, column=1)

if __name__ == "__main__":
    control_panel()
