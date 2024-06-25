import random 
import sys 
import time 
import os 
import pandas as pd 
import numpy as np 
import multiprocessing 
import matplotlib.pyplot as plt 
import py4j 
import nl4py 
from scipy.stats import chi2_contingency
from typing import List
import time

import skopt
from skopt import BayesSearchCV
from skopt import forest_minimize
from skopt import gp_minimize
from skopt.plots import plot_objective
print(f'skopt version: {skopt.__version__}')

import sklearn
from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
print(f'sklearn version: {sklearn.__version__}')


class simulation:
    def __init__(self, n_calls = 300, fname_NL_model = './Pulmonary_Histology_5_20_2024_Bill.nlogo'):
        self.fname_NL_model = fname_NL_model
        self.n_calls = n_calls
        self.start_time = time.time()
        nl4py.initialize("C:/Program Files/NetLogo 6.2.1")
        
    def est_ipf(self, parameterVal):

        # Loop through each value in the original list and duplicate it 24 times
        parameterVals = [parameterVal] * 24
        validation_distribution = [6,8,22] # distribution of fibrosis scores for people with IPF
        # define the model
        model = "./Pulmonary_Histology_5_20_2024_Bill.nlogo"

        # define simulation callback function that returns a list of setup commands
        def run_simulation(colDepDegTransVals) -> List[str]:    
            #Define parameters
            thy1ps_count = 90
            thy1ns_count = 10
            tnf_thresh = 50
            il1b_thresh = 50
            transition_threshold = colDepDegTransVals[2]
            tgfbint = 700
            coldep = colDepDegTransVals[0]
            coldeg = colDepDegTransVals[1]
            as_entry_threshold = 4    
            #Define setup commands
            setup_commands = ["reset-ticks", "setup", f"set thy1ps-count {thy1ps_count}", f"set thy1ns-count {thy1ns_count}",
            f"set set-tnf-thresh {tnf_thresh}", f"set set-il1b-thresh {il1b_thresh}",  
            f"set transition-threshold {transition_threshold}", f"set tgfbint {tgfbint}", 
            f"set coldep {coldep}", f"set coldeg {coldeg}", f"set as-entry-threshold {as_entry_threshold}"]
                
            return setup_commands    

        # metrics to report
        reporters = ["ticks", "fibrosis-score", "alv-count", "final-thy1n", "final-thy1p", "colcontent", "coldep", "coldeg", "TGFBint"]

        #use run_experiment to run n simulations and return results
        attempt = 0
        while attempt<5:
            try:
                results = nl4py.run_experiment(self.fname_NL_model, run_simulation, parameterVals, reporters, go_command = "repeat 365 [go]", stop_at_tick = 2)
                results_float = results["fibrosis-score"].to_numpy().astype(float)
                attempt = 10
            except:
                print(f'Attempt: {attempt}')
                attempt = attempt + 1
                
        if attempt==10:
            # if the simulation ran
            observed_distribution = [np.sum(results_float==0),np.sum(results_float==1),np.sum(results_float==2)]
            obs = np.array([observed_distribution, validation_distribution])
            stat = chi2_contingency(obs).statistic
        else:
            # if the simulation errored out 5 times
            print('WANRING: Simulation in NL4PY failed 5 times.')
            stat = 100
        elapsed_time = time.time()-self.start_time
        print(f'time={round(elapsed_time,2)}, statistic={round(stat,2)}, coldep={round(parameterVal[0],2)}, coldeg={round(parameterVal[1],2)}, transition_threshold = {round(parameterVal[2],2)}')
        return stat

    def sampling_minimization(self, method='RF'):

        bounds = [(0.0, 200.0), (0.0, 200.0), (0.0,0.2)]
        n_calls = self.n_calls

        if method=='RF':
            self.result = forest_minimize(
                self.est_ipf, bounds, n_calls=n_calls, base_estimator="ET", random_state=4
            )
        else:
            self.result = gp_minimize(
                self.est_ipf, bounds, n_calls=n_calls, random_state=4
            )
        
        return self.result
    
    def plot(self):
        _ = plot_objective(self, n_points=150)
    


