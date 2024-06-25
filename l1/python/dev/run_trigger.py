#!/usr/bin/env python
'''
date: 18 June 2024
author: jmgarriz

Script that runs the first implementation of the triggering algorithm for the full P-ONE detector.
As of 6/18/24 it only goes through "L0". 

'''
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units
from icecube.dataclasses import ModuleKey
from glob import glob
import argparse
import os
import numpy as np
from collections import defaultdict
from dataclasses import dataclass
import get_l0_trig as l0
import l1_lc_alg as l1

#add arguments for input, output, coincidence time window requirements and number of modules to be considered an event
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
parser.add_argument("-g", "--geo", default="/mnt/research/IceCube/PONE/jp_pone_sim/PONE_10String_7Cluster_standard_GCD.i3.gz",
                    help="I3 file containing the geometry of the simulation")
args = parser.parse_args()

#this function actually does all the icetray stuff 
def getData(frame):
    modules = defaultdict(list)
    if frame.Has("PMTResponse_nonoise"):
        pulsemap = frame['PMTResponse_nonoise']
        for omkey, p in pulsemap:
        #create a dictionary with (string, om) as keys and [pmt, pulses] as items
            key1 = ModuleKey(omkey[0], omkey[1])     
            modules[key1].append((omkey[2], p))
            #print(type(key1))
    #print(frame)
    #this might be inefficient and getting the geometry more times than necessary 
    if frame.Has("I3Geometry"):
        geometry = frame["I3Geometry"].omgeo
       #print(geometry)
        print("I got a geometry")
    
    #gather L0 triggers
    triggers = l0.findModuleMultiplicity(modules, args.window, args.moduleReq)
    #time order triggers
    ordered_triggers = l0.time_order(triggers)
    #for i in ordered_triggers:
        #print(i.multiplicity)
    test = l1.LC_reco_events(geometry, ordered_triggers)
    #print(test)

  
    

t = I3Tray()

t.AddModule("I3Reader", FilenameList=[args.geo, args.infile])
t.AddModule(getData, "getData", Streams = [icetray.I3Frame.DAQ])


t.Execute(6)
