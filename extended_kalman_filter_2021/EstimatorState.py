import numpy as np

"""
Class of all the 3 estimated states 

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
class EstimatorState:
    def __init__(self):
        self.Pm = np.empty((7,7))

        self.xm = np.empty(7)

        self.tm = np.empty(1)