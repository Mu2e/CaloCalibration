""" Mu2E Calorimeter Track event display version 1.3
    by Giacinto Boccia
    2024-08-13
    """

import ROOT as R
import numpy as np
#import concurrent.futures
from array import array
import crystalpos
#crystalpos.py is just a file with 2 np arrays with crystal x and y and the measure of the crystal side
CRYSTALS = 674

class Hit:
    def __init__(self, number, x, y, time, qval) -> None:
        self.n = number
        self.x = x
        self.y = y
        self.t = time
        self.q = qval

class Crystal:
    side = crystalpos.crys_side
    #Each crystal Should contain 2 hits

    def __init__(self, number : int, x : float, y : float) -> None:
        self.hit_arr = []
        self.n = number
        self.x = x
        self.y = y

    def angles(self):
        x = np.zeros(2)
        y = np.zeros(2)
        x[0] = self.x + self.side / 2
        y[0] = self.y + self.side / 2
        x[1] = self.x - self.side / 2
        y[1] = self.y - self.side / 2
        return x, y
    
    def force_new_hit(self, new_hit : Hit) -> int:
        #Returns the index where the hit is stored
        self.hit_arr.append(new_hit)
        return len(self.hit_arr) - 1
  
    def get_q(self) -> np.double:
        #If the crystal has multiple hits, an average is retunred
        if len(self.hit_arr):
            q_average = np.mean([hit.q for hit in self.hit_arr])
            return q_average
        else:
            return np.nan

    def test_new_hit(self, new_hit : Hit) -> bool:
        #Check if hit is inside the crystal, get angles
        x, y = self.angles()
        x_min = np.min(x)
        x_max = np.max(x)
        y_min = np.min(y)
        y_max = np.max(y)
        if x_min <= new_hit.x <= x_max and y_min <= new_hit.y <= y_max:
            self.force_new_hit(new_hit)
            return True
        else:
            return False
    
class Event:
    def __init__(self, event_number, tree_slice) -> None:
        self.crys_arr = []        
        self.fit_arr = []
        self.ev_num = event_number
        n_hits = tree_slice.nHits

        #prepare Crystals
        for i in range(CRYSTALS):
            self.crys_arr.append(Crystal(i, crystalpos.crys_x[i], crystalpos.crys_y[i]))

        #Load Hits
        x_arr = np.frombuffer(tree_slice.Xval)
        y_arr = np.frombuffer(tree_slice.Yval)
        t_arr = np.frombuffer(tree_slice.Tval)
        q_arr = np.frombuffer(tree_slice.Qval)
        for i in range(n_hits):
            #Loop over hits
            curr_hit = Hit(i, x_arr[i], y_arr[i], t_arr[i], q_arr[i])
            for crystal in self.crys_arr:
                #Loop over crystals
                if crystal.test_new_hit(curr_hit):
                    #When one of the crystals accepts the hit, go to next hit
                    break

    def event_prepare_hist(self) -> None:
        #This method prepares an histogram with the crystal q values
        self.crys_hist = R.TH2Poly()
        for crystal in self.crys_arr:
            boud_x, boud_y = crystal.angles()
            bin = self.crys_hist.AddBin(boud_x[0], boud_y[0], boud_x[1], boud_y[1])
            mean_q = crystal.get_q()
            if not np.isnan(mean_q):
                self.crys_hist.SetBinContent(bin, crystal.get_q())
    
    def event_fit(self) -> 'Event.Event_fit':
        #All stuff related to fitting events should be in the nested class Event_fit, this meas that you can fit the same event more than once (maybe with different parameters) without duplicating the event, just call this method mutiple times. All the fits created are available in a list at Event.fit_arr
        new_fit = Event.Event_fit(self)
        self.fit_arr.append(new_fit)
        return new_fit
    
    def event_draw(self, fits : bool = True ) -> None:
        #This method creates a new TCanvas for the whole event, and draws eveithing (hit, hinstogram, fits...)over it, if you want a TCanvas for a single fit, look for Event.Event_Fit.draw()
        ev_name = "Event " + str(self.ev_num)
        self.canvas = R.TCanvas(ev_name, ev_name, 1000, 1000)
        #R.gPad.SetGrid(1, 1)
        self.__histo_draw(ev_name)
        #This method can be invoked with no fits done, then it will draw just the detector activation, or with multiple fits, and it will overlay them
        if len(self.fit_arr) > 0 and fits:
            for fit in self.fit_arr:
                fit._Event_fit__fit_draw(ev_name)
        
        #Circes
        inner_c = R.TEllipse(0, 0, 374)
        inner_c.SetLineColor(2)
        inner_c.SetLineWidth(3)
        inner_c.SetFillStyle(0)
        inner_c.Draw('pl same')
        outer_c = R.TEllipse(0, 0, 660)
        outer_c.SetLineColor(2)
        outer_c.SetLineWidth(3)
        outer_c.SetFillStyle(0)
        outer_c.Draw('pl same')
        self.canvas.Draw()        

    def __histo_draw(self, ev_name : str) -> None:
        #Please use event_draw!
        if hasattr(self, 'crys_hist'):
            self.crys_hist.SetStats(0)
            self.crys_hist.SetTitle(ev_name)
            self.crys_hist.Draw('apl')  #'apl'
            self.crys_hist.Draw('zcol Cont0 same')
        else:
            print ("The event you are trying to draw has no crystal histogram! Try event_prepare_hist() before.")

    def centroid(self) -> tuple[np.double, np.double]:
        #Returns the event centroid
        x_arr = np.zeros(CRYSTALS)
        y_arr = np.zeros(CRYSTALS)
        q_arr = np.zeros(CRYSTALS)
        for i in range(CRYSTALS):
            x_arr[i] = self.crys_arr[i].x
            y_arr[i] = self.crys_arr[i].y
            q = self.crys_arr[i].get_q()
            if not np.isnan(q):
                q_arr[i] = q
        return np.average(x_arr, weights = q_arr), np.average(y_arr, weights = q_arr)

    class Event_fit:
        #Once you decide to fit an event, you instantiate this nested class, if new fit methods are developed, they should be placed in this class. To use them, call event_fit() on the Event, it returns an Event_fit object, if you vant to fit the same event vith different parametrs call event_fit() more than once
        def __init__(self, event : 'Event') -> None:
            self.event = event
            self.vertical = False

        def linear_fit(self, treshold : float) -> int:
            #The treshold applies on Qvals (mean value if crystal has 2 hits) and all hits over it get fitted linearly, returns the number of crystals considered
            n_sel = 0
            x_arr = array('f')
            y_arr = array('f')

            #Apply treshold
            for crystal in self.event.crys_arr:
                if crystal.get_q() > treshold:
                    x_arr.append(crystal.x)
                    y_arr.append(crystal.y)
                    n_sel +=1

            if n_sel > 1:
                #Put points in a TGraph
                self.n_sel = n_sel
                err_arr = array ('f', np.full(n_sel, crystalpos.crys_side / 2))               
                self.graph = R.TGraphErrors(n_sel, x_arr, y_arr, err_arr, err_arr)
                
                #Check if the points are on a verticl line
                if not(self.__is_vertical(x_arr, y_arr)):
                    self.fit = R.TF1("Linear", "pol1")
                    self.graph.Fit("Linear")
                    #Set the color here, so that different fits (differnet methods) can have different colors
                    self.fit.SetLineColor(4)
                    self.fit.SetLineWidth(2)
                else:
                    center = self.event.centroid()
                    self.fit = R.TLine(center[0], 660, center[0], -660)
                    self.fit.SetLineColor(4)
                    self.fit.SetLineWidth(2)
            return n_sel

        def __is_vertical(self, x_arr : array, y_arr : array) -> bool:
            #This method is called when fitting (it is suposed to be shared between differnet fit tecniques), if you want to access its result just read Event.Evnet_fit.vertical
            min_x = np.min(x_arr)
            max_x = np.max(x_arr)
            dx = max_x - min_x
            #If all xs are within a crystal side of the oters, the track is considerend vertical
            self.vertical = dx < crystalpos.crys_side
            return self.vertical

        def draw(self) -> None:
            #This method produces a TCanvas with the detector activation and the current fit (the one Event.Event_fit object on wich you call the method), if you are looking at a TCanvas with all the fits onverlayed over the calorimeter activation, or if you don't want to draw fits at all look at Event.event_draw()
            fit_name = "Event " + str(self.event.ev_num)
            self.canvas = R.TCanvas(fit_name, fit_name, 1000, 1000)
            #R.gPad.SetGrid(1, 1)
            self.event._Event__histo_draw(fit_name)
            self.__fit_draw(fit_name)
            
            #Circes
            inner_c = R.TEllipse(0, 0, 374)
            inner_c.SetLineColor(2)
            inner_c.SetLineWidth(3)
            inner_c.SetFillStyle(0)
            inner_c.Draw('pl same')
            outer_c = R.TEllipse(0, 0, 660)
            outer_c.SetLineColor(2)
            outer_c.SetLineWidth(3)
            outer_c.SetFillStyle(0)
            outer_c.Draw('pl same')
            self.canvas.Draw()
                      
        def __fit_draw(self, name: str) -> None:
            #Please use draw() instead!
            if hasattr(self, 'graph'):
                #Style
                self.graph.SetMarkerColor(1)
                self.graph.SetMarkerSize(2)
                self.graph.GetXaxis().SetTitle('X (mm)')
                self.graph.GetYaxis().SetTitle('Y (mm)')
                self.graph.GetXaxis().SetLimits(-660, 660)
                self.graph.GetYaxis().SetRangeUser(-650, 650)
                self.graph.Draw('* same')
                if hasattr(self, 'fit'):
                    self.fit.Draw('same')
            else:
                print("You are drawing a fit that does not exist! Try linear_fit() [or other] before.")
            
def tree_loader(tree, n_max = 0):
    #Loads a tree (up to n_max, if set) in an list of Event objects
    ev_arr = []
    n_events = 0
    for slice in tree:
        event = Event(n_events, slice)
        ev_arr.append(event)
        n_events += 1
        #Just a break if you want to run on a fraction of the tree
        if n_max != 0 and  n_events >= n_max:
            break
    return ev_arr


if __name__ == '__main__':
    TRESHOLD = 4000
    BINS = 500
    EVLIMIT = 200

    #Open tree
    file_name = input("File to open:")
    file = R.TFile.Open(file_name)
    tree = file.sidet

    #Load events
    ev_arr = tree_loader(tree, EVLIMIT)

    while True:
        ev_num = int(input("Event to display:"))
        event = ev_arr[ev_num]
        event.event_prepare_hist()
        fit = event.event_fit()
        n_points = fit.linear_fit(TRESHOLD)
        print ("Hits above treshold:", n_points)
        fit.draw()