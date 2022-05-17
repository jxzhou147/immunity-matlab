from SALib.analyze import sobol
from SALib.analyze import sobol
import numpy as np

from sobol_sample import problem

# load Q values for txt files
Q_V = np.loadtxt("QV.txt", float)
Q_M = np.loadtxt("QM.txt", float)
Q_K = np.loadtxt("QK.txt", float)

# compute sensitivity indices
Si_V = sobol.analyze(problem, Q_V)
Si_M = sobol.analyze(problem, Q_M)
Si_K = sobol.analyze(problem, Q_K)
