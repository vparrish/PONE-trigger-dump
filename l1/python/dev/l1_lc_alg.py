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
#this LC algorithm now matches LC alg we've been referencing 
def light_cone(dist):
    #define constant variables 
    c = 0.299792458 #m/ns
    n = 1.34 #index of refraction
    d_atten = 80.5 #m, not sure what this exactly is
    theta_c = math.acos(1/n) #radians, not really sure what value this is
   

    #lower limit
    
    small_dist_lim = 2*d_atten*math.sin(theta_c) 
    if dist < small_dist_lim:
        tmin = 0
    elif dist >= small_dist_lim:
        tmin_int = (1/c)*np.sqrt(dist**2 - small_dist_lim**2)
        #print(tmin_int)
        tmin_large = (d_atten/c)*((1/n) + np.sqrt((dist/d_atten)**2 - math.sin(theta_c)**2) -n)
        #print(tmin_large)
        if tmin_int < tmin_large:
            tmin = tmin_int
            #print("intermediate")
        else:
            tmin = tmin_large
            #print("large")
    #upper limit
            
    tmax = (n*dist)/c
    if dist <= d_atten:
        tmax =(n*dist)/c
    else:
        tmax = (d_atten/c)*(np.sqrt((dist/d_atten)**2 - math.sin(theta_c)**2) - (1/n) +n)

    pos_max =  tmax
    pos_min = tmin
    neg_max = -tmin
    neg_min = -tmax

    #returns 4 values for the pos/neg max/min values of the windows
    return pos_max, pos_min, neg_max, neg_min
    
def distance(x_t, y_t, z_t, x_m, y_m, z_m):
    mag = np.sqrt((x_t-x_m)**2 + (y_t-y_m)**2 + (z_t-z_m)**2)
    return mag


def LC_reco_events(geometry, triggers, max_dist):
    #insert something to go thru all triggers and then find distance between the trigger pulse and the 
    #every optical module
    d_max = max_dist
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
                    pos_max, pos_min, neg_max, neg_min = light_cone(dist)
                    
                    #add something in here to then from those time windows and the module distances, see if there is a pulse we want to look at ? 
                    pos_tw = [t.time+pos_min, t.time+pos_max]
                    neg_tw = [t.time - neg_min, t.time - neg_max]
                    

                    #include the omkey of the neighbor module so that we can then determine what omkeys and time windows to search from the pmtresponse no noise frame
                    light_cone_data.append((pos_max, pos_min, neg_max, neg_min, dist, t.time, omkey))
                    
            temp = ModuleKey(omkey[0], omkey[1])
    light_cone_df = pd.DataFrame(light_cone_data, columns=['pos_max', 'pos_min', 'neg_max', 'neg_min', 'distance', 'time', 'omkey'])
    return light_cone_df

def reduce_windows(light_cone_sampled):
    #this function will go thru and select which time windows for which omkeys we should search so we can avoid duplicate searches 
    pos_max = light_cone_sampled['pos_max']
    pos_min = light_cone_sampled['pos_min']
    neg_max = light_cone_sampled['neg_max']
    neg_min = light_cone_sampled['neg_min']
    dist = light_cone_sampled['distance']
    time = light_cone_sampled['time']
    omkey = light_cone_sampled['omkey']

    
    #for i in range(len(omkey)):
    

def plotlc(light_cone_df):
    # Sample 1% of the data
    light_cone_sampled = light_cone_df.sample(frac=1)

    pos_max = light_cone_sampled['pos_max']
    pos_min = light_cone_sampled['pos_min']
    neg_max = light_cone_sampled['neg_max']
    neg_min = light_cone_sampled['neg_min']
    dist = light_cone_sampled['distance']
    time = light_cone_sampled['time']

    # Creating subplots
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))  # 2 Rows, 1 Column

    axs[0].scatter(dist, pos_max, marker='.', alpha=0.5, label='positive max', color = "teal")
    axs[0].scatter(dist, pos_min, marker='.', alpha=0.5, label='positive min', color  = "deeppink")
    axs[0].set_title('Time Window Positive')
    axs[0].set_xlabel('Module separation')
    axs[0].set_ylabel('Positive time delay')
    #axs[0].set_xscale('log')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].scatter(dist, neg_max, marker='.', alpha=0.5, label='negative max', color= "teal")
    axs[1].scatter(dist, neg_min, marker='.', alpha=0.5, label='negative min', color = "deeppink")
    axs[1].set_title('Time Window Negative')
    axs[1].set_xlabel('Module separation')
    axs[1].set_ylabel('Negative time delay')
    #axs[1].set_xscale('log')
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()  # Adjust layout to not overlap
    plt.savefig("lc_alg_plt.png")


def eventClustering(neighbors):
    potential_events= []
    #perform operation to combine time windows and modules. this might be its own function this will be an annoying function to write :) 
    #group these together as events
    return potential_events

#next add a function to go back into the simulation and look to see if we see a pulse in those time windows 

    
