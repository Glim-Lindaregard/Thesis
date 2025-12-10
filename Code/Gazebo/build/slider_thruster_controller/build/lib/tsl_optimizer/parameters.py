import numpy as np 


MAX_FORCE = 0.7
MAX_TAU = 0.0098

nt = 8              # Number of thrusters
ni = 3              # Number of inputs
A = np.array([     # Thruster matrix
    [-1,     0,     1,     0,     1,     0,    -1,     0   ],
    [ 0,     1,     0,     1,     0,    -1,     0,    -1   ],
    [-0.14,  0.14,  0.14, -0.14, -0.14,  0.14,  0.14, -0.14],
], dtype=np.float64)

P = np.identity(8)
q = np.zeros((8, 1))
G = np.identity(8)
h = np.array([
        [MAX_FORCE], [MAX_FORCE], [MAX_FORCE], [MAX_FORCE], [MAX_FORCE], [MAX_FORCE], [MAX_FORCE], [MAX_FORCE]], 
        dtype=np.float64) 