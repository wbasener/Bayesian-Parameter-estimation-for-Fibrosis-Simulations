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


from skopt import BayesSearchCV
from skopt import forest_minimize
from skopt import gp_minimize
from skopt.plots import plot_objective

from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split

plt.rcParams["figure.figsize"] = (12,12)

nl4py.initialize("C:/Program Files/NetLogo 6.2.1")



def est_ipf(parameterVal):

    # Loop through each value in the original list and duplicate it 12 times
    parameterVals = [parameterVal] * 36
    validation_distribution = [6,8,22] # distribution of fibrosis scores for people with IPF
    # define the model
    model = "./Pulmonary_Histology_5_20_2024_Bill.nlogo"

    #define simulation callback function that returns a list of setup commands
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

    #metrics to report
    reporters = ["ticks", "fibrosis-score", "alv-count", "final-thy1n", "final-thy1p", "colcontent", "coldep", "coldeg", "TGFBint"]

    #use run_experiment to run n simulations and return results
    check = False
    while check==False:
        try:
            results = nl4py.run_experiment(model, run_simulation, parameterVals, reporters, go_command = "repeat 365 [go]", stop_at_tick = 3)
            results_float = results["fibrosis-score"].to_numpy().astype(float)
            check = True
        except:
            check = False
    observed_distribution = [np.sum(results_float==0),np.sum(results_float==1),np.sum(results_float==2)]
    obs = np.array([observed_distribution, validation_distribution])
    stat = chi2_contingency(obs).statistic
    
    print(f'statistic={stat}, coldep={parameterVal[0]}, coldeg={parameterVal[1]}, transition_threshold = {parameterVal[2]}')
    return stat


bounds = [(0.0, 200.0), (0.0, 200.0), (0.0,0.2)]
n_calls = 300


result_RG = forest_minimize(
    est_ipf, bounds, n_calls=n_calls, base_estimator="ET", random_state=4
)
_ = plot_objective(result, n_points=150)
plt.title('X_0 = ColDep, X_1 = ColDeg, X_2 = transition_threshold')
plt.savefig('partial_dpendence_RF')


result_GP = gp_minimize(
    est_ipf, bounds, n_calls=n_calls, random_state=4
)
_ = plot_objective(result, n_points=150)
plt.title('X_0 = ColDep, X_1 = ColDeg, X_2 = transition_threshold')
plt.savefig('partial_dpendence_GP')

