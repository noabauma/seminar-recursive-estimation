import numpy as np
import multiprocessing
import time

np.set_printoptions(precision=4)

from run import run_estimator
from Simulator import Simulator
from SimulationConst import SimulationConst
from EstimatorConst import EstimatorConst



def main():
    trackErrorNorm = run_estimator(doplot=True, seed=113)

    print('trackErrorNorm: ', trackErrorNorm)


if __name__ == "__main__":
    main()