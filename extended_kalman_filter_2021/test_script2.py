import numpy as np
import sys
import time

"""
#to check which modules are loaded
modulenames = set(sys.modules) & set(globals())
allmodules = [sys.modules[name] for name in modulenames]
print(allmodules)
"""

from EstimatorConst import EstimatorConst
from SimulationConst import SimulationConst
from Simulator import *
from run import *

np.set_printoptions(precision=4, threshold=sys.maxsize)

def main():

    N = 1000
    
    time_arr = np.zeros(N)
    trackErrorNorm_arr = np.zeros(N)
    angularErrorNorm_arr = np.zeros(N)
    velocityErrorNorm_arr = np.zeros(N)
    windErrorNorm_arr = np.zeros(N)
    biasErrorNorm_arr = np.zeros(N)

    for i in range(N):
        start = time.time()
        trackErrorNorm,angularErrorNorm,velocityErrorNorm,windErrorNorm,biasErrorNorm,initialEst,initialVar = run(doplot=False, seed=i+42)
        
        time_arr[i]              = time.time() - start
        trackErrorNorm_arr[i]    = trackErrorNorm
        angularErrorNorm_arr[i]  = angularErrorNorm
        velocityErrorNorm_arr[i] = velocityErrorNorm
        windErrorNorm_arr[i]     = windErrorNorm
        biasErrorNorm_arr[i]     = biasErrorNorm

    print('trackErrorNorm:')
    print(np.mean(trackErrorNorm_arr))
    print(np.var(trackErrorNorm_arr, ddof=1))

    print('angularErrorNorm:')
    print(np.mean(angularErrorNorm_arr))
    print(np.var(angularErrorNorm_arr, ddof=1))

    print('velocityErrorNorm:')
    print(np.mean(velocityErrorNorm_arr))
    print(np.var(velocityErrorNorm_arr, ddof=1))

    print('windErrorNorm:')
    print(np.mean(windErrorNorm_arr))
    print(np.var(windErrorNorm_arr, ddof=1))

    print('biasErrorNorm:')
    print(np.mean(biasErrorNorm_arr))
    print(np.var(biasErrorNorm_arr, ddof=1))

    print('time:')
    print(np.mean(time_arr))
    print(np.var(time_arr, ddof=1))
    
    
    

if __name__ == "__main__":
    main()