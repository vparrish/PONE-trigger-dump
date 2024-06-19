'''
date: 19 June 2024
author: jmgarriz

Script to perform LC algorithm calculation event clustering
!!WORK IN PROGRESS, INCOMPLETE!!
'''

from glob import glob
import argparse
import os
import numpy as np
from collections import defaultdict
from dataclasses import dataclass
import math

from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units


parser = argparse.ArgumentParser(
    description="collects all the l0 triggers and outputs to some format not yet determined :)")

parser.add_argument("-i", "--infile", default="/mnt/research/IceCube/PONE/jp_pone_sim/pmtsim/GenerateSingleMuons_39_pmtsim.i3.zst",
                    help="input .gz files")
parser.add_argument("-o", "--outfile", default="./out.i3",
                    help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-w", "--window", default=10,#ns
                    help="length of coincidence time window")
parser.add_argument("-m", "--moduleReq", default=2,
                    help="Minimum number of modules which must have PEs for an event to be considered")

args = parser.parse_args()

@dataclass
class LC_event:
    def __init__(self, module, posWindow, negWindow, trigModule):
        self.module = module
        self.posWindow = posWindow
        self.negWindow = negWindow
        self.trigModule = trigModule


#initial start on implementing light cone algorithm function
#made to be independent of icetray stuff and just take already calculated distance and time
#untested atm
def light_cone(dist, t_initial):
    #define constant variables 
    c = 0.299792458 #m/ns
    n = 1.34 #index of refraction
    d_atten = 25 #m, not sure what this exactly is 
    theta_c = 40.5 #degrees?, not really sure what value this is
    d_max = 500 #m, we need some maximum distance at which we just won't even bother calculating

    #lower limit
    if dist > d_max:
        return 0
    else:
        small_dist_lim = 2*d_atten*math.sin(math.radians(theta_c)) 
        if dist < small_dist_lim:
            tmin = 0
        elif dist >= small_dist_lim:
            tmin_int = (1/c)*np.sqrt(dist**2 - small_dist_lim**2)
            tmin_large = (d_atten/c)*((1/n) + np.sqrt((dist/d_atten)**2 - math.sin(math.radians(theta_c))**2) -n)
        if tmin_int < tmin_large:
            tmin = tmin_int
        else:
            tmin = tmin_large
    #upper limit
            
    tmax = (n*dist)/c
    if dist <= (tmax*c)/n:
        tmax =(n*dist)/c
    elif dist < d_max:
        tmax = (d_atten/c)*np.sqrt((dist/d_atten)**2 - math.sin(math.radians(theta_c))**2 - (1/n) +n)

    #positive and negative time windows
    
    tw_p = [t_initial + tmin, t_initial +tmax]
    tw_n = [t_initial -  tmin, t_initial - tmax]

    time_windows = [tw_p, tw_n]
    return time_windows

def dist(x_t, y_t, z_t, x_m, y_m, z_m):
    mag = np.sqrt((x_t-x_m)**2 + (y_t-y_m)**2 + (z_t-z_m)**2)
    return mag


def LC_reco_events(frame, triggers):
    #insert something to go thru all triggers and then find distance between the trigger pulse and the 
    #every optical module 
    key = frame["I3OMGeoMap"]
    neighbors = []
    for i in triggers:
        for j in key:
            t = i.module
            print(t)
            #distance = dist(t.pos.x, t.pos.y, t.pos.z, j.pos.x, j.pos.y, j.pos.z)
            #print(distance)
            #compute distance here- if same module skip 
            #perform LC alg based on distance and time of the trigger pulse 
            #if LC function does not result in return zero, add to LC_event 
            #append LC_events to a neighbors list per moduletrigger 
            #print('placeholder so vscode stops giving me issues') 
        
    return neighbors 
        
def eventClustering(neighbors):
    potential_events= []
    #perform operation to combine time windows and modules. this might be its own function this will be an annoying function to write :) 
    #group these together as events
    return potential_events

#next add a function to go back into the simulation and look to see if we see a pulse in those time windows 

    
