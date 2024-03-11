from glob import glob
import pandas as pd
from icecube import dataio, phys_services
from icecube import dataclasses
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="performs cascade reconstruction on simulated events using spline table")

parser.add_argument("-i", "--infile", default="./test_input.i3",
                    help="input .gz files")


def get_time_residuals():
    frames = []
    infiles = sorted(glob(args.infile))
    for file in infiles:
        print(file)
        infile = dataio.I3File(file)
        while infile.more():
            try:
                frames.append(infile.pop_DAQ())
            except:
                continue
        infile.close()
    print(frames)