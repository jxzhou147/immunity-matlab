from SALib.sample import saltelli
import numpy as np
import pandas as pd
import math

par_name = [
    'a_MI', 'a_KI', 'gamma', 'c_M_IL6', 'c_M_IL10', 'c_M_CCL2', 'c_I_CXCL5'
]

par_IAV = pd.read_csv('par_IAV.txt', header=None)
par_IAV = par_IAV[1].to_numpy()
par_consider = par_IAV[[4, 6, 1, 15, 19, 20, 23]]

ratio = 2

par_bounds = []
for i in range(len(par_consider)):
    par_bounds.append([
        par_consider[i] - math.log10(ratio),
        par_consider[i] + math.log10(ratio)
    ])

# Define the model inputs
problem = {'num_vars': 7, 'names': par_name, 'bounds': par_bounds}

# Generate samples
par_values = saltelli.sample(problem, 128)
np.savetxt("par_values.txt", par_values)