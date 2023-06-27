import numpy as np
import multiprocessing
import time

#np.random.seed(1)
#np.set_printoptions(precision=4)

from run import run_estimator
from Simulator import Simulator
from SimulationConst import SimulationConst
from EstimatorConst import EstimatorConst

def worker(x):
    start = time.time()
    trackErrorNorm = run_estimator(doplot=False, seed=int(x+1))

    return trackErrorNorm, time.time() - start, x



def main():

    pool = multiprocessing.Pool(processes=12)
    pool_data = pool.map(worker, range(100))

    pool_data = np.array(pool_data)

    trackErrorNorms = pool_data[:,0]
    Exctimes = pool_data[:,1]

    print(pool_data)

    print('trackErrorNorms:')
    print(np.mean(trackErrorNorms))
    print(np.var(trackErrorNorms, ddof=1))

    print('Execution times:')
    print(np.mean(Exctimes))
    print(np.var(Exctimes, ddof=1))



    

if __name__ == "__main__":
    main()