from SALib.sample import saltelli
import numpy as np

par_name = [
    'a_MI', 'a_KI', 'gamma', 'c_M_IL6', 'c_M_IL10', 'c_M_CCL2', 'c_I_CXCL5'
]

par_bounds = [[-6.2309, -1.5577], [-3.5506, -0.8876], [1.6010, 6.4039],
              [-3.1257, -0.7814], [-1.9605, -0.4901], [-1.3263, -0.3316],
              [-0.3835, -0.0959]]

# Define the model inputs
problem = {'num_vars': 7, 'names': par_name, 'bounds': par_bounds}

# Generate samples
par_values = saltelli.sample(problem, 64)
np.savetxt("par_values.txt", par_values)