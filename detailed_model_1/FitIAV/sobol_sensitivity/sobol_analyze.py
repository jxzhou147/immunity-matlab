from SALib.analyze import sobol
from SALib.analyze import sobol
import numpy as np

from sobol_sample import problem

# load Q values from txt files
# Q_V = np.loadtxt("QV.txt", float)
# Q_M = np.loadtxt("QM.txt", float)
# Q_K = np.loadtxt("QK.txt", float)

# load y values from txt files
y_V = np.loadtxt("y_V.txt", float)
y_M = np.loadtxt("y_M.txt", float)
y_K = np.loadtxt("y_K.txt", float)

# compute sensitivity indices
# Si_V = sobol.analyze(problem, Q_V, print_to_console=True)
# Si_M = sobol.analyze(problem, Q_M)
# Si_K = sobol.analyze(problem, Q_K)

Si_V = [sobol.analyze(problem, y) for y in y_V.T]
Si_M = [sobol.analyze(problem, y) for y in y_M.T]
Si_K = [sobol.analyze(problem, y) for y in y_K.T]

S1s_V = np.array([s['S1'] for s in Si_V])
S1s_M = np.array([s['S1'] for s in Si_M])
S1s_K = np.array([s['S1'] for s in Si_K])

STs_V = np.array([s['ST'] for s in Si_V])
STs_M = np.array([s['ST'] for s in Si_M])
STs_K = np.array([s['ST'] for s in Si_K])

# write sobol indices to files
np.savetxt("S1s_V.txt", S1s_V)
np.savetxt("S1s_M.txt", S1s_M)
np.savetxt("S1s_K.txt", S1s_K)

np.savetxt("STs_V.txt", STs_V)
np.savetxt("STs_M.txt", STs_M)
np.savetxt("STs_K.txt", STs_K)