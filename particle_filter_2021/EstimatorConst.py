import numpy as np

""" 
 Define the physical constants that are available to the estimator.


 Class:
 Recursive Estimation
 Spring 2021
 Programming Exercise 2

 --
 ETH Zurich
 Institute for Dynamic Systems and Control
 Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
 hofermat@ethz.ch
 csferrazza@ethz.ch

 --
 Revision history
 [21.04.19, CS]    first version by Matthias & Carlo
 [29.03.21, TB]    second version by Thomas
"""
class EstimatorConst:
    def __init__(self):
        """ The wall contour - (x,y) coordinates of corners p1 to p10 as in Table 1 """
        self.contour = np.array([[0.50, 0.00],
                            [2.50, 0.00],
                            [2.50, 1.50],
                            [3.00, 2.00],
                            [2.00, 3.00],
                            [1.25, 2.25],
                            [1.00, 2.50],
                            [0.00, 2.50],
                            [0.00, 0.50],
                            [0.50, 0.50]])
                    
        """ Initialization """
        self.pA = [1.1, 0.6] # Center point pA of the initial position distribution
        self.pB = [1.8, 2.0] # Center point pB of the initial position distribution
        self.d = 0.2  # Radius of shaded regions for initialization

        self.phi_0 = np.pi/4.0 # Initial heading is uniformly distr. in [-phi_0,phi_0]

        self.l = 0.2  # Uniform distribution parameter of points p9 and p10

        """ Noise properties """
        # process noise
        self.sigma_phi = 0.05 # Parameter for process noise v_phi
        self.sigma_f = 0.01 # Parameter for process noise v_f

        # measurement noise
        self.epsilon = 0.01 # Parameter for measurement noise w
