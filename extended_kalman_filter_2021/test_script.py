import numpy as np
import sys

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
    estConst = EstimatorConst()
    simConst = SimulationConst()

    N = 1000
    
    """
    trackErrorNorm_arr = np.zeros(N)
    angularErrorNorm_arr = np.zeros(N)
    velocityErrorNorm_arr = np.zeros(N)
    windErrorNorm_arr = np.zeros(N)
    biasErrorNorm_arr = np.zeros(N)

    for i in range(N):
        trackErrorNorm,angularErrorNorm,velocityErrorNorm,windErrorNorm,biasErrorNorm,initialEst,initialVar = run(doplot=False, seed=i+42)
        
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
    """
    
    tm_means = np.zeros(N)
    tm_vars  = np.zeros(N)

    state_means = np.zeros((N,5))
    state_vars = np.zeros((N,5))

    wind_means = np.zeros(N)
    wind_vars  = np.zeros(N)

    drift_means = np.zeros(N)
    drift_vars  = np.zeros(N)

    input_means = np.zeros((N,2))
    input_vars  = np.zeros((N,2))

    sense_means = np.zeros((N,5))
    sense_vars  = np.zeros((N,5))


    for i in range(N):
        tm, state, wind, drift, input, sense = Simulator(simConst)

        tm_means[i] = np.mean(tm)
        tm_vars[i]  = np.var(tm, ddof=1)

        state_means[i,:] = np.mean(state, axis=0)
        state_vars[i,:] = np.var(state, axis=0, ddof=1)

        wind_means[i] = np.mean(wind)
        wind_vars[i]  = np.var(wind, ddof=1)

        drift_means[i] = np.mean(drift)
        drift_vars[i]  = np.var(drift, ddof=1)

        input_means[i,:] = np.mean(input, axis=0)
        input_vars[i,:]  = np.var(input, axis=0, ddof=1)

        sense[sense == np.inf] = np.nan
        sense[sense == -np.inf] = np.nan

        sense_means[i,:] = np.nanmean(sense, axis=0)
        sense_vars[i,:]  = np.nanvar(sense, axis=0, ddof=1)

    print('tm_means:')
    print(np.mean(tm_means))
    print(np.var(tm_means, ddof=1))
    print('tm_vars:')
    print(np.mean(tm_vars))
    print(np.var(tm_vars, ddof=1))

    print('state_means:')
    print(np.mean(state_means, axis=0))
    print(np.var(state_means, axis=0, ddof=1))
    print('state_vars:')
    print(np.mean(state_vars, axis=0))
    print(np.var(state_vars, axis=0, ddof=1))

    print('wind_means:')
    print(np.mean(wind_means))
    print(np.var(wind_means, ddof=1))
    print('wind_vars:')
    print(np.mean(wind_vars))
    print(np.var(wind_vars, ddof=1))

    print('drift_means:')
    print(np.mean(drift_means))
    print(np.var(drift_means, ddof=1))
    print('drift_vars:')
    print(np.mean(drift_vars))
    print(np.var(drift_vars, ddof=1))

    print('input_means:')
    print(np.mean(input_means, axis=0))
    print(np.var(input_means, axis=0, ddof=1))
    print('input_vars:')
    print(np.mean(input_vars, axis=0))
    print(np.var(input_vars, axis=0, ddof=1))

    print('sense_means:')
    print(np.mean(sense_means, axis=0))
    print(np.var(sense_means, axis=0, ddof=1))
    print('sense_vars:')
    print(np.mean(sense_vars, axis=0))
    print(np.var(sense_vars, axis=0, ddof=1))

    """
    tm, state, wind, drift, input, sense = Simulator(simConst)

    print("\ntm")
    print(tm.shape)
    print(np.mean(tm))
    print(np.var(tm))

    print("\nstate")
    print(state.shape)
    print(np.mean(state, axis=0))
    print(np.var(state, axis=0))

    print("\nwind")
    print(wind.shape)
    print(np.mean(wind))
    print(np.var(wind))

    print("\ndrift")
    print(drift.shape)
    print(np.mean(drift))
    print(np.var(drift))

    print("\ninput")
    print(input.shape)
    print(np.mean(input, axis=0))
    print(np.var(input, axis=0))

    sense[sense == np.inf] = np.nan
    sense[sense == -np.inf] = np.nan
    print("\nsense")
    print(sense.shape)
    #print(sense)
    print(np.nanmean(sense, axis=0))
    print(np.nanvar(sense, axis=0))
    """

if __name__ == "__main__":
    main()