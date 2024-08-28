import l1_lc_alg as l1
import get_l0_trig as l0
import numpy as np
import math 
import matplotlib.pyplot as plt
from icecube.dataclasses import ModuleKey
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="collects all the l0 triggers and outputs to some format not yet determined :)")
parser.add_argument("-g", "--geo", default="/mnt/research/IceCube/PONE/jp_pone_sim/PONE_10String_7Cluster_standard_GCD.i3.gz",
                    help="I3 file containing the geometry of the simulation")
args = parser.parse_args()


def show_geo(frame):
    if frame.Has("I3Geometry"):
        geometry = frame["I3Geometry"].omgeo
        #print(type(geometry))
        print("I got a geometry")
    #triggers is output as a list of module triggers which includes the module key, multiplicity, and the interaction time
    #need to generate a list of module triggers that we can put into the geometry and perform the LC calculation on it
    #times for these should all be zero (abs time doesn't matter)
    #mult doesnt matter either just put 2 in 
    mks = []
    x = []
    y = []
    temp = ModuleKey(999,999)
    for omkey, pos in geometry:
        mk = ModuleKey(omkey[0], omkey[1])
        #print(mk)
        if mk !=temp:
            string = mk[0]
            mks.append(string)
            x.append(pos.position[0])
            y.append(pos.position[1])
        temp = mk
    print(len(x))
    print(len(mks))
    #MT = l0.ModuleTrigger(mk, 2, 0) 
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    for i, txt in enumerate(mks):
        if i%20 ==0:
            ax.text(x[i], y[i], txt)
    plt.savefig("geo.png")
    
#based on function from internet
def merge_overlap(arr):
    # Sort intervals based on start values
    arr.sort(key=lambda x: x[0])

    res_idx = 0  # Index of the last merged interval

    for i in range(1, len(arr)):
      
        # If current interval overlaps with 
        # the last merged interval
        if arr[res_idx][1] >= arr[i][0]:
            arr[res_idx][1] = max(arr[res_idx][1], arr[i][1])
            #print("merging")
        else:
            # Move to the next interval
            #print("moving on")
            res_idx += 1
            arr[res_idx] = arr[i]

    # Return the size of the merged intervals
    return res_idx + 1


def remove_duplicates(geometry,max_dist, rate):
    #sample trigger 
    triggers = []
    t = l0.ModuleTrigger(ModuleKey(1, 1), 3, 0)
    triggers.append(t)
    data_frame = l1.LC_reco_events(geometry, triggers, max_dist)
    pos_max = data_frame['pos_max']
    pos_min = data_frame['pos_min']
    neg_max = data_frame['neg_max']
    neg_min = data_frame['neg_min']
    pos_arr = []
    neg_arr = []
    for i in range(len(pos_max)):
        pos_arr.append([pos_min[i], pos_max[i]])
        neg_arr.append([neg_min[i], neg_max[i]])
    print(len(pos_arr))
    new_size_pos = merge_overlap(pos_arr)
    for i in range(new_size_pos):
        print(f"[{pos_arr[i][0]}, {pos_arr[i][1]}]", end=" ")
    print(new_size_pos)

    new_size_neg = merge_overlap(neg_arr)
    for i in range(new_size_neg):
        print(f"[{neg_arr[i][0]}, {neg_arr[i][1]}]", end=" ")



def run_study(frame):
    if frame.Has("I3Geometry"):
        geometry = frame["I3Geometry"].omgeo
    remove_duplicates(geometry, 100, 2000)


t = I3Tray()

t.AddModule("I3Reader", FilenameList=[args.geo])
t.AddModule(run_study, "run_study", Streams = [icetray.I3Frame.Geometry])


t.Execute(6)