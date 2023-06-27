import time
import numpy as np
import matplotlib.pyplot as plt
import copy

from SimulationConst import SimulationConst
from EstimatorConst import EstimatorConst
from EstimatorState import EstimatorState

from Simulator import Simulator
from Estimator import Estimator

#Helper functions
two_pi = 2.0*np.pi
def wrapTo2Pi(x):
    tmp = x % two_pi
    positiveInput = x > 0.0
    tmp[(tmp == 0.0) & positiveInput] = two_pi #useless feature in my opinion
    return tmp

#returns row-wise dot product
def multi_dot(x,y):
    if x.shape != y.shape:
        raise Exception('x & y have to be the same shape! x.shape = ' + str(x.shape) +  'y.shape = ' + str(y.shape))

    dim = x.shape[0]
    tmp = np.empty(dim)
    for i in range(dim):
        tmp[i] = x[i,:].dot(y[i,:])

    return tmp

"""
 Main function for Extended Kalman Filter programming exercise.
 It is used to simulate the true model, call the estimator and show the 
 results.


 Class:
 Recursive Estimation
 Spring 2022
 Programming Exercise 1

 --
 ETH Zurich
 Institute for Dynamic Systems and Control
 Raffaello D'Andrea, Enrico Mion, Bhavya Sukhija, Jin Cheng
 enmion@ethz.ch
 sukhijab@ethz.ch
 jicheng@student.ethz.ch

"""
def run(simConst=SimulationConst(), estConst=EstimatorConst(), doplot=True, seed=0):
    """ Setup """

    # Set the random number generator state.
    # Set != to == or give a seed not equal to 0 
    # to make results reproducable. This setting was used to generate
    # the plot in the problem description.
    if seed != 0:
        np.random.seed(seed)

    # TODO: maybe check for the amount of libraries used?
    # like: numpy, scipy, casadi, matplotlib
    # do it here!

    """ Simulation """
    # The function 'Simulator' simulates the boat dynamics and generates
    # measurements.
    tm, state, wind, drift, input, sense = Simulator(copy.deepcopy(simConst))

    # state contains rolling out of [p_x, p_y, s_x, s_y, phi]: time index from 1 to N

    # input contains rolling out of [u_t, u_r]: time index from 1 to N-1

    # sense contains rolling out of [z_a, z_b, z_c, z_g, z_n] (contains inf):
    # time index from 1 to N, with first measurement inf (to indicate no
    # measurement at initial step)

    """ Run the Estimator """
    N = simConst.N

    # Initialize the estimator.  
    estState  = EstimatorState()
    posEst    = np.zeros((N,2))
    linVelEst = np.zeros((N,2))
    oriEst    = np.zeros(N)
    windEst   = np.zeros(N)
    driftEst  = np.zeros(N)
    posVar    = np.zeros((N,2))
    linVelVar = np.zeros((N,2))
    oriVar    = np.zeros(N)
    windVar   = np.zeros(N)
    driftVar  = np.zeros(N)

    # Call the estimator for each time step.

    # Initialization
    posEst[0,:],linVelEst[0,:],oriEst[0],windEst[0],driftEst[0],posVar[0,:],linVelVar[0,:],oriVar[0],windVar[0],driftVar[0],estState = Estimator(estState, [], [], tm[0], copy.deepcopy(estConst))
    
    # Remaining steps
    tstart = time.time()
    for n in range(1,N):
        posEst[n,:],linVelEst[n,:],oriEst[n],windEst[n],driftEst[n],posVar[n,:],linVelVar[n,:],oriVar[n],windVar[n],driftVar[n],estState = Estimator(estState, input[n-1,:], sense[n,:], tm[n], copy.deepcopy(estConst))
    
    tend = time.time()
    tEstAvg = (tend - tstart)/(N-1)

    if doplot:
        print('Done. Average estimator single update computation time: ', tEstAvg)

    """ The results """
    # Calculate the total tracking error as an indicator of your estimator's performance.
    trackError = np.hstack((state[:,0] - posEst[:,0], state[:,1] - posEst[:,1]))
    trackErrorNorm = np.sqrt(trackError.dot(trackError)/N)

    unitVectorTrue = np.vstack((np.cos(state[:,4]), np.sin(state[:,4]))).T
    unitVectorEst = np.vstack((np.cos(oriEst),np.sin(oriEst))).T
    angularError = wrapTo2Pi(np.arccos(multi_dot(unitVectorTrue,unitVectorEst)))
    angularErrorNorm = np.sqrt(angularError.dot(angularError)/N)

    velocityError = np.hstack((state[:,2] - linVelEst[:,0], state[:,3] - linVelEst[:,1]))
    velocityErrorNorm = np.sqrt(velocityError.dot(velocityError)/N)

    unitVectorTrue = np.array([np.cos(wind),np.sin(wind)]).T
    unitVectorEst = np.array([np.cos(windEst),np.sin(windEst)]).T
    windError = wrapTo2Pi(np.arccos(multi_dot(unitVectorTrue,unitVectorEst)))
    windErrorNorm = np.sqrt(windError.dot(windError)/N)

    unitVectorTrue = np.array([np.cos(drift),np.sin(drift)]).T
    unitVectorEst = np.array([np.cos(driftEst),np.sin(driftEst)]).T
    biasError = wrapTo2Pi(np.arccos(multi_dot(unitVectorTrue,unitVectorEst)))
    #biasError =  wrapTo2Pi(drift) - wrapTo2Pi(driftEst)
    biasErrorNorm = np.sqrt(biasError.dot(biasError)/N)

    initialEst = [posEst[0,:],linVelEst[0,:],oriEst[0],driftEst[0]]
    initialVar = [posVar[0,:],linVelVar[0,:],oriVar[0],driftVar[0]]

    # Plots of the results.
    if doplot:
        fig, axs = plt.subplots(3,3)

        i = 0
        row = int(i/3.0)
        col = i % 3

        # planar position plot
        axs[row,col].axis('equal')
        axs[row,col].plot(state[:,0], state[:,1],'r.', label='true state')
        axs[row,col].plot(posEst[:,0], posEst[:,1],'b.', label='estimate')
        axs[row,col].set_xlabel('x position [m]')
        axs[row,col].set_ylabel('y position [m]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # x position over time
        axs[row,col].plot(tm, state[:,0],'r.', label='true state')
        axs[row,col].plot(tm, posEst[:,0],'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('x position [m]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # y position over time
        axs[row,col].plot(tm, state[:,1],'r.', label='true state')
        axs[row,col].plot(tm, posEst[:,1],'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('y position [m]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
            
        # orientatin over time
        axs[row,col].plot(tm, state[:,4],'r.', label='true state')
        axs[row,col].plot(tm, oriEst,'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('orientation [rad]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # x velocity over time
        axs[row,col].plot(tm, state[:,2],'r.', label='true state')
        axs[row,col].plot(tm, linVelEst[:,0],'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('x velocity  [m/s]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # y velocity over time
        axs[row,col].plot(tm, state[:,3],'r.', label='true state')
        axs[row,col].plot(tm, linVelEst[:,1],'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('y velocity [m/s]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # drift over time
        axs[row,col].plot(tm, drift,'r.', label='true state')
        axs[row,col].plot(tm, driftEst,'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('gyro drift [rad]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3
        
        # wind speed over time
        axs[row,col].plot(tm, wind, 'r.', label='true state')
        axs[row,col].plot(tm, windEst,'b.', label='estimate')
        axs[row,col].set_xlabel('time [s]')
        axs[row,col].set_ylabel('wind direction [rad]')
        axs[row,col].legend()
        axs[row,col].grid()

        i += 1
        row = int(i/3.0)
        col = i % 3

        axs[row,col].set_visible(False)

        #plt.savefig("plots.png")
        plt.show()
    
    return trackErrorNorm,angularErrorNorm,velocityErrorNorm,windErrorNorm,biasErrorNorm,initialEst,initialVar