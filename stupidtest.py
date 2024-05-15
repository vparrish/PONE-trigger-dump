from glob import glob
import pandas as pd
from icecube import dataio, phys_services
from icecube import dataclasses
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="performs cascade reconstruction on simulated events using spline table")

parser.add_argument("-i", "--infile", default="GenerateSingleMuons_39_pmtsim.i3.zst")
args = parser.parse_args()

def get_time_residuals():
    frames = []
    infiles = sorted(glob(args.infile))
    for file in infiles:
        print(file)
        infile = dataio.I3File(file)
        while infile.more():
            try:
                frames.append(infile.pop_daq())
                print("yay")
            except:
                continue
        infile.close()
    print(frames)

get_time_residuals()




'''
modules = defaultdict(lambda: defaultdict(PulseQueue))
        for p_key, p_value in pulsemap:
            module_key = (p_key[0], p_key[1])
            pmt_pulses = {p_key[2]: PulseQueue(p_value[0], p_value[-1])}
            modules[module_key].update(pmt_pulses)
        print(modules)

'''