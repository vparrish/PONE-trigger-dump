import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import prelim_study as ps
from icecube.dataclasses import ModuleKey
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units
import argparse

df = pd.read_pickle("calls_df.pkl")
#print(df)

calls = df['call_num']
string = df['string_num']
om = df['om_num']

#print(calls)

plt.hist(calls, bins = 50)
plt.title("average number of calls (10kHz noise rate)")
plt.ylabel("number of calls")
plt.savefig("calls_hist.png")
print(np.average(calls))

parser = argparse.ArgumentParser(
    description="collects all the l0 triggers and outputs to some format not yet determined :)")
parser.add_argument("-g", "--geo", default="/mnt/research/IceCube/PONE/jp_pone_sim/PONE_10String_7Cluster_standard_GCD.i3.gz",
                    help="I3 file containing the geometry of the simulation")
args = parser.parse_args()

