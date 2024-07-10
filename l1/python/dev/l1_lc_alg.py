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
import matplotlib.pyplot as plt
import pandas as pd
from icecube.dataclasses import ModuleKey
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

class trig_geo:
    def __init__(self, mk, geo, time):
        self.mk = mk
        self.geo = geo
        self.time = time
                 

#initial start on implementing light cone algorithm function
#made to be independent of icetray stuff and just take already calculated distance and time
#untested atm
def light_cone(dist, t_initial, d_max):
    #define constant variables 
    c = 0.299792458 #m/ns
    n = 1.34 #index of refraction
    d_atten = 25 #m, not sure what this exactly is 
    theta_c = 40.5 #degrees?, not really sure what value this is
   

    #lower limit
    
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
        
    tw_p_bounds = [t_initial + tmin, t_initial +tmax]
    tw_n_bounds = [t_initial -  tmin, t_initial - tmax]

    tw_p_central = (tw_p_bounds[0] + tw_p_bounds[1])/2
    tw_n_central = (tw_n_bounds[0] + tw_n_bounds[1])/2
    time_windows = [tw_p_bounds, tw_n_bounds]

    return tw_p_central, tw_n_central, tw_p_bounds, tw_n_bounds

def distance(x_t, y_t, z_t, x_m, y_m, z_m):
    mag = np.sqrt((x_t-x_m)**2 + (y_t-y_m)**2 + (z_t-z_m)**2)
    return mag


def LC_reco_events(geometry, triggers):
    #insert something to go thru all triggers and then find distance between the trigger pulse and the 
    #every optical module
    d_max = 500
    neighbors = []
    trigger_geometries = []
    for t in triggers:
        p = t.module
        for omkey, pos in geometry:
            mk = ModuleKey(omkey[0], omkey[1])
            if mk == p:
                t_geo = trig_geo(mk, pos, t.time)
                trigger_geometries.append(t_geo)
    #this temp thing is here so that we don't repeat the distance calculation for every pmt on same module
    temp = ModuleKey(999,999)
    light_cone_data = []
    for t in trigger_geometries:
        for omkey, pos in geometry:
            mk = ModuleKey(omkey[0], omkey[1])
            if mk == t.mk:
                #print("got the same modules")
                continue
            elif mk != temp:
                #print(temp.om)
                t_x = t.geo.position.x
                t_y = t.geo.position.y
                t_z = t.geo.position.z
                m_x = pos.position.x
                m_y = pos.position.y
                m_z = pos.position.z
                dist = distance(t_x, t_y, t_z, m_x, m_y, m_z)
                #print(dist)
                if dist <= d_max:
                    tw_p_central, tw_n_central, tw_p_bounds, tw_n_bounds = light_cone(dist, t.time, d_max)
                    #so this is printing some kinds of time windows now....
                    #unsure if they are correct or not and what format we want them in tbh
                    # print(lc)
                    light_cone_data.append((tw_p_central, tw_n_central, tw_p_bounds, tw_n_bounds, dist, t.time))
                    
            temp = ModuleKey(omkey[0], omkey[1])

    
            #perform LC alg based on distance and time of the trigger pulse 
            #if LC function does not result in return zero, add to LC_event 
            #append LC_events to a neighbors list per moduletrigger 
            #print('placeholder so vscode stops giving me issues') 
    light_cone_df = pd.DataFrame(light_cone_data, columns=['tw_p_central', 'tw_n_central', 'tw_p_bounds', 'tw_n_bounds', 'distance', 'time'])
    return neighbors, light_cone_df

def plotlc(light_cone_df):
    light_cone_df = light_cone_df.sample(frac=0.1, random_state=1)
    # Extracting central values and bounds for tw_p
    tw_p_central = light_cone_df['tw_p_central']
    # Setting bounds as error
    tw_p_bounds = np.array(light_cone_df['tw_p_bounds'].tolist())
    tw_p_errors = [abs(tw_p_bounds[:, 0]), abs(tw_p_bounds[:, 1])]

    # Extracting central values and bounds for tw_n
    tw_n_central = light_cone_df['tw_n_central']
    tw_n_bounds = np.array(light_cone_df['tw_n_bounds'].tolist())
    tw_n_errors = [abs(tw_n_bounds[:, 0]), abs(tw_n_bounds[:, 1])]

    # Creating subplots
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))  # 2 Rows, 1 Column

    # Plotting tw_p vs distance with error bars
    axs[0].errorbar(light_cone_df['distance'], tw_p_central, yerr=tw_p_errors, fmt='o', alpha=0.5, label='central tw pos', ecolor='lightgray', elinewidth=3, capsize=0)
    axs[0].set_title('Time Window Positive')
    axs[0].set_xlabel('Module separation')
    axs[0].set_ylabel('Positive time delay')
    axs[0].set_xscale('log')
    axs[0].grid(True)
    axs[0].legend()

    # Plotting tw_n vs distance with error bars
    axs[1].errorbar(light_cone_df['distance'], tw_n_central, yerr=tw_n_errors, fmt='o', alpha=0.5, color='red', label='central tw neg', ecolor='lightgray', elinewidth=3, capsize=0)
    axs[1].set_title('Time Window Negative')
    axs[1].set_xlabel('Module separation')
    axs[1].set_ylabel('Negative time delay')
    axs[1].set_xscale('log')
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()  # Adjust layout to not overlap
    plt.show()


def eventClustering(neighbors):
    potential_events= []
    #perform operation to combine time windows and modules. this might be its own function this will be an annoying function to write :) 
    #group these together as events
    return potential_events

#next add a function to go back into the simulation and look to see if we see a pulse in those time windows 

    
