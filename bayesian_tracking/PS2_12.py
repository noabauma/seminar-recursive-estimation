import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

np.random.seed(1313)  #random fixed seed (for reproducibility)

"""
NOTE: This is a python version of the original code which runs matlab.
To run the code: 
install using: "pip3 install numpy matplotlib" 
and then run "python3 PS2_12.py"
"""
def main():
    """ Configuration Constants """

    # number of simulation steps
    T = 100

    # number of discrete steps around circle
    N = 100

    # actual probability of going CCW
    PROB = 0.55

    # model of probability of ging CCW
    PROB_MODEL = PROB
    #PROB_MODEL = 0.45

    # Location of distance sensor, as a multiple of the circle radius.  Can be
    # less than 1 (inside circle), but must be positive (WLOG).
    SENSE_LOC = 2.0

    # The sensor error is modeled as additive (a time of flight sensor, for
    # example), uniformly distributed around the actual distance.  The units
    # are in circle radii.  
    ERR_SENSE = 0.50

    # Model of what the sensor error is
    ERR_SENSE_MODEL = ERR_SENSE
    #ERR_SENSE_MODEL = 0.45

    """ Initialization """

    # W(k,i) denotes the probability that the object is at location i at time
    # k, given all measurements up to, and including, time k.  At time 0, this
    # is initialized to 1/N, all positions are equally likely.
    W = np.zeros((T+1,N))
    W[0,:] = 1.0/N

    # The intermediate prediction weights, initialize here for completeness.
    # We don't keep track of their time history.
    predictW = np.zeros(N)

    # The initial location of the object, an integer between 0 and N-1.
    loc = np.zeros(T+1)
    loc[0] = round(N/4)

    """ Simulation """
    for t in range(1, T+1):
        """ Simulate System """

        # Process dynamics. With probability PROB we move CCW, otherwise CW
        if np.random.rand() < PROB:
            loc[t] = (loc[t-1] + 1) % N
        else:
            loc[t] = (loc[t-1] - 1) % N

        # The physical location of the object is on the unit circle
        xLoc = np.cos(2*np.pi * loc[t]/N)
        yLoc = np.sin(2*np.pi * loc[t]/N)

        # Can calculate the actual distance to the object
        dist = np.sqrt((SENSE_LOC - xLoc)**2 +  yLoc**2)

        # Corrupt the distance by noise
        dist += ERR_SENSE * 2.0 * (np.random.rand() - 0.5)

        """ Update Estimator """
        
        
        for i in range(N):
            # Prediction Step.  Here we form the intermediate weights which capture
            # the pdf at the current time, but not using the latest measurement.
            predictW[i] = PROB_MODEL*W[t-1, (i-1)%N] + (1.0-PROB_MODEL)*W[t-1, (i+1)%N]

            # Fuse prediction and measurement.  We simply scale the prediction step
            # weight by the conditional probability of the observed measurement
            # at that state.  We then normalize.
            xLocHypo = np.cos(2*np.pi*i/N)
            yLocHypo = np.sin(2*np.pi*i/N)

            distHypo = np.sqrt((SENSE_LOC - xLocHypo)**2 + yLocHypo**2)

            if np.abs(dist - distHypo) < ERR_SENSE_MODEL:
                condProb = 1.0/(2.0*ERR_SENSE_MODEL)
            else:
                condProb = 0.0

            W[t,i] = condProb * predictW[i]

        # Normalize the weights.  If the normalization is zero, it means that
        # we received an inconsistent measurement.  We can either use the old
        # valid data, re-initialize our estimator, or crash. To be as
        # robust as possible, we simply re-initialize the estimator.
        normConst = np.sum(W[t,:])
        
        # Comment this line if we want to allow the program to crash.
        if normConst > 1e-6:
            W[t,:] /= normConst
        else:
            W[t,:] = W[0,:]

    """  Visualize the results """
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    xVec = np.linspace(0, (N-1)/N, N)
    yVec = np.arange(T+1)
    X, Y = np.meshgrid(xVec, yVec)

    ax.set_xlabel('POSITION x(k)/N')
    ax.set_ylabel('TIME STEP k')

    ax.plot_surface(X, Y, W, cmap=cm.coolwarm)
    ax.plot(loc/N, yVec, np.ones(T+1)*np.max(W))

    plt.show()  # Windows users probably have to use: plt.savefig("plot.png")


if __name__ == "__main__":
    main()