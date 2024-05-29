import random 
import sys 
import time 
import os 
import pandas as pd 
import numpy as np 
import multiprocessing 
import matplotlib.pyplot as plt 
#from deap import base, creator, tools, algorithms 
import py4j 
import nl4py 
from scipy.stats import chi2_contingency
from typing import List
import site
print(site.getsitepackages())

#initialize netlogo
print("Initializing NL4Py...")
nl4py.initialize("C:/Program Files/NetLogo 6.2.1")

# Baseline: (From Rikard et al)
colDepBase = 0.0056
colDegBase = 1 # From Rikard et al: 0.0035

#sample parameter values for the experiment
scl = 2000 # generate values up to scl times the bas value
nuSamples = 4 # number of samples to generate
depDeg_ratio = 200 # max ratio of colDep to colDeg, that is, colDeg < colDep < depDeg_ratio*colDeg
colDepDegVals = []
for i in range(nuSamples):
    colDeg = np.round(np.random.uniform(colDegBase, scl*colDegBase), 3) # colDegBase < colDeg < scl*colDegBase
    colDep = np.round(np.random.uniform(colDeg, depDeg_ratio*colDeg), 3) # colDeg < colDep < depDeg_ratio*colDeg
    colDepDegVals.append([colDep, colDeg])

# Loop through each value in the original list and duplicate it 36 times
colDepDegVals_dup = []
for value in colDepDegVals:
    colDepDegVals_dup.extend([value] * 36)

#check
print("Parameter values sampled:")
for v in colDepDegVals:
    print(f'[colDep, colDeg] = {v}')

#define the model
model = "./Pulmonary_Histology_5_20_2024_Bill.nlogo"

#define simulation callback function that returns a list of setup commands
def run_simulation(colDepDegVals) -> List[str]:    
    #Define parameters
    thy1ps_count = 90
    thy1ns_count = 10
    tnf_thresh = 50
    il1b_thresh = 50
    transition_threshold = 0.1 #false is equal to 0
    tgfbint = 700
    coldep = colDepDegVals[0]
    coldeg = colDepDegVals[1]
    as_entry_threshold = 4    
    #Define setup commands
    setup_commands = ["reset-ticks", "setup", f"set thy1ps-count {thy1ps_count}", f"set thy1ns-count {thy1ns_count}",
     f"set set-tnf-thresh {tnf_thresh}", f"set set-il1b-thresh {il1b_thresh}",  
     f"set transition-threshold {transition_threshold}", f"set tgfbint {tgfbint}", 
     f"set coldep {coldep}", f"set coldeg {coldeg}", f"set as-entry-threshold {as_entry_threshold}"]
        
    return setup_commands    

#metrics to report
reporters = ["ticks", "fibrosis-score", "alv-count", "final-thy1n", "final-thy1p", "colcontent", "coldep", "coldeg", "TGFBint"]

#use run_experiment to run n simulations and return results
#this configuration will run 3 sets of 365 steps and report metrics at 365, 730, and 1095
#stop_at_tick=365 will output a FS for every single tick point
print("Running simulations and collecting results...")
results = nl4py.run_experiment(model, run_simulation, colDepDegVals_dup, reporters, go_command = "repeat 365 [go]", stop_at_tick = 3)
print("Finished running simulations...")

#save the results in a data frame
print("Saving results...")
print(results)
df = pd.DataFrame(results)
print(df)
df.to_csv("results.csv", index=False)




