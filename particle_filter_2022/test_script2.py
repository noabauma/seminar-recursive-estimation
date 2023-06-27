import numpy as np
import multiprocessing
import time

#np.set_printoptions(precision=4)

from run import run_estimator
from Simulator import Simulator
from SimulationConst import SimulationConst
from EstimatorConst import EstimatorConst



def main():
    N = 1000

    trackError_arr = []
    time_arr = []
    for i in range(N):
        tstart = time.time()
        tmp1 = run_estimator(doplot=False, seed=2*i + 4369)
        tmp2 = time.time() - tstart

        trackError_arr.append(tmp1)
        time_arr.append(tmp2)

        if i % 10 == 0:
            print('i = ', i, flush=True)

            print('trackErrorNorm:')
            print(np.mean(trackError_arr), flush=True)
            print(np.var(trackError_arr, ddof=1), flush=True)

            print('execution time:')
            print(np.mean(time_arr), flush=True)
            print(np.var(time_arr, ddof=1), flush=True)


    print("final")
    print('trackErrorNorm:')
    print(np.mean(trackError_arr))
    print(np.var(trackError_arr, ddof=1))

    print('execution time:')
    print(np.mean(time_arr))
    print(np.var(time_arr, ddof=1))


if __name__ == "__main__":
    main()