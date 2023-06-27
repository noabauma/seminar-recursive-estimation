import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import copy

from SimulationConst import SimulationConst
from EstimatorConst import EstimatorConst

from Simulator import Simulator
from Estimator import Estimator

from Animation import Animation

# trackErrorNorm=run_estimator(SimulationConstVar1(),EstimatorConstVar1(),doplot,seed)
#
# Main function for the Particle Filter programming exercise.
# It is used to simulate the true model, call the estimator and show the 
# results.
#
#
# Class:
# Recursive Estimation
# Spring 2021
# Programming Exercise 2
#
# --
# ETH Zurich
# Institute for Dynamic Systems and Control
# Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
# hofermat@ethz.ch
# csferrazza@ethz.ch
def run_estimator(simConst=SimulationConst(), estConst=EstimatorConst(), doplot=True, seed=0):
    """ Setup """
    # Set the random number generator state.
    # Uncomment to make results reproducable.
    if seed != 0:
        np.random.seed(seed)

    # TODO: maybe check for the amount of libraries used?
    # like: numpy, scipy, casadi, matplotlib
    # do it here!

    """ Simulation """
    # The function 'Simulator' simulates the robot dynamics and generates
    # measurements.
    km, state, input, sense = Simulator(copy.deepcopy(simConst))

    # state contains rolling out of [x_r, y_r, phi, kappa] for k={1,N}
    # input contains rolling out of [u_f, u_phi] for k={1,N-1}
    # sense contains rolling out of [z] for k={1,N}, where first measurement is
    # infinity (to indicate no measurement at initial step)

    """ Run the Estimator """
    N = simConst.N 

    # Initialize the estimator.
    posterior = Estimator([], [], [], copy.deepcopy(estConst), km[0])

    # Get number of particles:
    N_particles = posterior.x_r.shape[0]

    # initialize storage arrays for x, y, phi and kappa:
    x_r = np.zeros((N, N_particles))
    x_r[0,:] = posterior.x_r
    y_r = np.zeros((N, N_particles))
    y_r[0,:] = posterior.y_r
    phi = np.zeros((N, N_particles))
    phi[0,:] = posterior.phi
    kappa = np.zeros((N, N_particles))
    kappa[0,:] = posterior.kappa

    # Call the estimator for the remaining time steps and measure the total 
    # execution time (including the updates of the storage arrays):
    if doplot:
        print('Generating data...')
    
    nextDispFracComplete = 0.1
    tstart = time.time()
    for k in range(1,N):
        if doplot:
            if ( k/N ) > nextDispFracComplete:
                nextDispFracComplete += 0.1
                print("Progress {:2.0%}".format(k/N), end="\r")
        
        posterior = Estimator(copy.deepcopy(posterior), sense[k,:], input[k-1,:], copy.deepcopy(estConst), km[k])
        # update storage arrays:
        x_r[k,:]   = posterior.x_r
        y_r[k,:]   = posterior.y_r
        phi[k,:]   = posterior.phi
        kappa[k,:] = posterior.kappa
        
    tEstAvg = (time.time() - tstart)/(N-1)
    
    if doplot:
        print('Done. Average estimator single update computation time: ', tEstAvg)
        print('Making plots.')

    """ Plotting """
    if doplot:
        ani = Animation(N_particles, x_r, y_r, phi, kappa, state, sense, simConst, fps=120)
        
        #ani.save('test_ani.mp4', fps=120)
        
        plt.show()

    """ The results """
    # initialize tracking error
    trackError = np.zeros(N)

    # Calculate the total tracking error as an indicator of your estimator's performance.
    for k in range(N):
        trackError[k] = np.sqrt(np.mean((x_r[k,:] - state[k,0])**2 + (y_r[k,:] - state[k,1])**2))
    
    return np.sqrt(trackError.dot(trackError)/N)