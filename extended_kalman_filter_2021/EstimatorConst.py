import numpy as np

"""
Class of all the physical constants that are available to the estimator

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
class EstimatorConst:
    def __init__(self):
        """ Boat dynamics constants """
        self.dragCoefficientHydr = 0.1 #C_d,h
        self.dragCoefficientAir = 0.06 #C_d,a
        self.rudderCoefficient = 2.0 #C_r
        self.windVel = 0.75 # C_w

        """ Radio measuremtn constants """
        self.pos_radioA = [-1000, 1000] #[x_a, y_a]
        self.pos_radioB = [2000, 0] #[x_b, y_b]
        self.pos_radioC = [0, 2000] #[x_c, y_c]

        """ Noise properties """
        # process noise
        # defined by its variance Q_
        self.DragNoise = 0.1 # Q_d
        self.RudderNoise = 0.01 # Q_r
        self.WindAngleNoise = 0.01 # Q_rho
        self.GyroDriftNoise = 0.01 # Q_b

        # measurement noise
        # normally distributed with zero mean and variance \sigma_
        self.DistNoiseA = 20.0 # DistNoiseA = \sigma_a^2
        self.DistNoiseB = 20.0 # DistNoiseB = \sigma_b^2
        self.DistNoiseC = 5.0 # DistNoiseC = \sigma_c^2

        self.GyroNoise = 0.01 # \sigma_g^2

        self.CompassNoise = 0.5 #\sigma_n^2

        """ Starting point """
        # The boat start with equal probility in a circle of radius R0 around the 
        # origin
        self.StartRadiusBound = 10.0 # R_0

        # The initial orientation is uniformly distributed in the range
        # [-\bar{\phi}, \bar{\phi}] 
        self.RotationStartBound = np.pi/8.0 # \bar{\phi}

        # The initial wind direction is uniformly distributed in the range
        # [-\bar{\rho}, \bar{\rho}] 
        self.WindAngleStartBound = np.pi/16.0 # \bar{\rho}

        # The initial gyro drift is exactly 0
        self.GyroDriftStartBound = 0.0
