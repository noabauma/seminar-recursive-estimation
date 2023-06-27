import numpy as np
from scipy.integrate import solve_ivp #ode45 imitation

from EstimatorConst import EstimatorConst 
from EstimatorState import EstimatorState

"""
The estimator.

The function is initialized for tm == 0, otherwise the estimator does an
iteration step (compute estimates for the time step k).

Inputs:
   estState        previous estimator state (time step k-1)
                   May be defined by the user (for example as a struct).
   actuate         control input u(k-1), [1x2]-vector
                   actuate[0]: u_t, thrust command
                   actuate[1]: u_r, rudder command
   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
                   measurement
                   sense[0]: z_a, distance measurement a
                   sense[1]: z_b, distance measurement b
                   sense[2]: z_c, distance measurement c
                   sense[3]: z_g, gyro measurement
                   sense[4]: z_n, compass measurement
   tm              time t_k, scalar
                   If tm==0 initialization, otherwise estimator
                   iteration step.
   estConst        estimator constants (as in EstimatorConst.m)

Outputs:
   posEst          position estimate (time step k), [1x2]-vector
                   posEst[0]: p_x position estimate
                   posEst[1]: p_y position estimate
   linVelEst       velocity estimate (time step k), [1x2]-vector
                   linVelEst[0]: s_x velocity estimate
                   linVelEst[1]: s_y velocity estimate
   oriEst          orientation estimate (time step k), scalar
   windEst         wind direction estimate (time step k), scalar
   driftEst        estimate of the gyro drift b (time step k), scalar
   posVar          variance of position estimate (time step k), [1x2]-vector
                   posVar[0]: x position variance
                   posVar[1]: y position variance
   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
                   linVelVar[0]: x velocity variance
                   linVelVar[1]: y velocity variance
   oriVar          variance of orientation estimate (time step k), scalar
   windVar         variance of wind direction estimate(time step k), scalar
   driftVar        variance of gyro drift estimate (time step k), scalar
   estState        current estimator state (time step k)
                   Will be input to this function at the next call.


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

--
Revision history
[24.04.18, MH]    first version by Matthias & Carlo
[07.05.19, CS]    updated for 2019
[29.03.21, TB]    updated for 2021
[09.04.22, NB]    written into python & updated for 2022
"""
def Estimator(estState: EstimatorState, actuate, sense, tm, estConst: EstimatorConst):

    """ Initialization """
    if tm == 0:
        # Do the initialization of your estimator here!
        
        # Replace the following:
        posEst = [0,0] # initial pos means are zero
        linVelEst = [0,0] # initial vel means are zero
        oriEst = 0 # initial rot mean is zero
        windEst = 0 # initial wind direction is zero
        driftEst = 0 # initial drift is zero
        
        # initial position variance is calculated analytically
        posVar = [1/4*estConst.StartRadiusBound**2,1/4*estConst.StartRadiusBound**2]
        linVelVar = [0,0] # initial velocity variance is zero
        # initial orientation variance is calculated analytically
        oriVar = 1/3*estConst.RotationStartBound**2
        windVar = 1/3*estConst.WindAngleStartBound**2
        driftVar = 0 # initial drift variance is zero
        
        # estimator variance init (initial posterior variance)
        estState.Pm = np.diag(np.hstack((posVar, linVelVar, oriVar, driftVar, windVar)))
        # estimator state (pos, rot, drift)
        estState.xm = np.hstack((posEst, linVelEst, oriEst, driftEst, windEst))
        # time of last update
        estState.tm = tm

        return posEst,linVelEst,oriEst,windEst,driftEst,posVar,linVelVar,oriVar,windVar,driftVar,estState

    """ Estimator iteration. """
    # x = [p_x p_y s_x s_y phi b rho]' in R^7
    # v = [v_d, v_r, v_b v_rho]' in R^4
    # w = [w_a, w_b, w_c, w_g, w_n]' in R^5

    # get time np.since last estimator update
    dt = tm - estState.tm
    estState.tm = tm # update measurement update time

    # process noise Q
    Q = np.diag([estConst.DragNoise, estConst.RudderNoise, estConst.GyroDriftNoise, estConst.WindAngleNoise])

    """ prior update """

    # helper function
    def dynamics(t,x):    
        # helper vars:
        thrust = np.tanh(actuate[0])
        c_dh = estConst.dragCoefficientHydr
        c_da = estConst.dragCoefficientAir
        c_r  = estConst.rudderCoefficient
        c_w  = estConst.windVel
        
        # state ODE
        dx = np.zeros(56)
        dx[0] = x[2]                                     # dp_x/dt
        dx[1] = x[3]                                     # dp_y/dt
        dx[2] = np.cos(x[4])*(thrust-c_dh*(x[2]**2+x[3]**2))-c_da*(x[2]-c_w*np.cos(x[6]))*np.sqrt((x[2]-c_w*np.cos(x[6]))**2+(x[3]-c_w*np.sin(x[6]))**2)   # ds_x/dt
        dx[3] = np.sin(x[4])*(thrust-c_dh*(x[2]**2+x[3]**2))-c_da*(x[3]-c_w*np.sin(x[6]))*np.sqrt((x[2]-c_w*np.cos(x[6]))**2+(x[3]-c_w*np.sin(x[6]))**2)   # ds_y/dt
        dx[4] = c_r*actuate[1]                           # dphi/dt
        dx[5] = 0                                        # db/dt
        dx[6] = 0                                        # drho/dt
        
        # linearize for variance update
        A = np.array([[0, 0,                                                                                                                                                                                                       1,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0],
                      [0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       1,                                            0, 0,                                                                                                                                                                                                 0],
                      [0, 0, - c_da*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5 - 2*c_dh*x[2]*np.cos(x[4]) - (c_da*(x[2] - c_w*np.cos(x[6]))*(2*x[2] - 2*c_w*np.cos(x[6])))/(2*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5),                                                                - 2*c_dh*x[3]*np.cos(x[4]) - (c_da*(x[2] - c_w*np.cos(x[6]))*(2*x[3] - 2*c_w*np.sin(x[6])))/(2*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5), -np.sin(x[4])*(thrust - c_dh*(x[2]**2 + x[3]**2)), 0, (c_da*c_w*(x[2] - c_w*np.cos(x[6]))*(x[3]*np.cos(x[6]) - x[2]*np.sin(x[6])))/((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5 - c_da*c_w*np.sin(x[6])*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5],
                      [0, 0,                                                                - 2*c_dh*x[2]*np.sin(x[4]) - (c_da*(x[3] - c_w*np.sin(x[6]))*(2*x[2] - 2*c_w*np.cos(x[6])))/(2*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5), - c_da*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5 - 2*c_dh*x[3]*np.sin(x[4]) - (c_da*(x[3] - c_w*np.sin(x[6]))*(2*x[3] - 2*c_w*np.sin(x[6])))/(2*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5),  np.cos(x[4])*(thrust - c_dh*(x[2]**2 + x[3]**2)), 0, c_da*c_w*np.cos(x[6])*((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5 + (c_da*c_w*(x[3] - c_w*np.sin(x[6]))*(x[3]*np.cos(x[6]) - x[2]*np.sin(x[6])))/((x[2] - c_w*np.cos(x[6]))**2 + (x[3] - c_w*np.sin(x[6]))**2)**0.5],
                      [0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0],
                      [0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0],
                      [0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0]])

        L = np.array([[0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [-np.cos(x[4])*c_dh*(x[2]**2+x[3]**2), 0, 0, 0],
                      [-np.sin(x[4])*c_dh*(x[2]**2+x[3]**2), 0, 0, 0],
                      [0, dx[4], 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
        
        P = np.reshape(x[7:],(7,7))
        
        # variance ODE
        dP = np.matmul(A,P) + np.matmul(A,P.T) + np.matmul(np.matmul(L,Q),L.T)
        # combine
        dx[7:] = dP.T.flatten()

        return dx

    # use ODE to push mean and variance forward:
    y0 = np.hstack((estState.xm, estState.Pm.T.flatten()))

    sol = solve_ivp(dynamics, [0.0, dt], y0, t_eval=[0.0, dt])

    xsim = sol.y.transpose()

    xp = xsim[-1,:7]
    Pp = np.reshape(xsim[-1,7:],(7,7))

    """ measurement update """
    # sensor constants:
    x_a = estConst.pos_radioA[0]
    y_a = estConst.pos_radioA[1]
    x_b = estConst.pos_radioB[0]
    y_b = estConst.pos_radioB[1]
    x_c = estConst.pos_radioC[0]
    y_c = estConst.pos_radioC[1]

    # Build measurement update matrices (assume distance C meas. is available)
    H = np.array([[(xp[0]-x_a)/((xp[0]-x_a)**2+(xp[1]-y_a)**2)**0.5, 
                   (xp[1]-y_a)/((xp[0]-x_a)**2+(xp[1]-y_a)**2)**0.5, 0, 0, 0, 0, 0], # z_a
                  [(xp[0]-x_b)/((xp[0]-x_b)**2+(xp[1]-y_b)**2)**0.5,
                   (xp[1]-y_b)/((xp[0]-x_b)**2+(xp[1]-y_b)**2)**0.5, 0, 0, 0, 0, 0], # z_b
                  [(xp[0]-x_c)/((xp[0]-x_c)**2+(xp[1]-y_c)**2)**0.5,
                   (xp[1]-y_c)/((xp[0]-x_c)**2+(xp[1]-y_c)**2)**0.5, 0, 0, 0, 0, 0], # z_c
                  [0, 0, 0, 0, 1, 1, 0],  # z_g
                  [0, 0, 0, 0, 1, 0, 0]]) # z_n
    # M:
    M = np.array([[1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 1]])
    
    # R:
    R = np.diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise])
    # z:
    z = sense[:5]

    # Nonlinear measurement prediction
    h_k = np.array([((xp[0]-x_a)**2+(xp[1]-y_a)**2)**0.5,
                    ((xp[0]-x_b)**2+(xp[1]-y_b)**2)**0.5,
                    ((xp[0]-x_c)**2+(xp[1]-y_c)**2)**0.5,
                    xp[4]+xp[5],
                    xp[4]])

    # If distance A,B & C meas. is not available, reduce H, M, R, z, h_k matrices
    if np.isinf(sense[2]):
        H = H[[3,4],:]
        M = M[3:5,3:5]
        R = R[3:5,3:5]
        z = z[[3,4]]
        h_k = h_k[[3,4]]

    if z.size != 0: #TODO: how can z.size == 0?! try to understand
        # do update only if there was a measurement, otherwise, just keep 
        # prior update
        temp1 = np.matmul(Pp,H.T)
        temp2 = np.matmul(np.matmul(M,R),M.T)
        temp3 = np.matmul(np.matmul(H,Pp),H.T) + temp2
        K = np.matmul(temp1, np.linalg.inv(temp3))

        estState.xm = xp + np.matmul(K, (z - h_k).T)
        # Joseph form for EKF
        temp4 = np.identity(7) - np.matmul(K,H)
        estState.Pm = np.matmul(np.matmul(temp4, Pp), temp4.T) + np.matmul(np.matmul(K, temp2), K.T)
    else:
        # if there was no meas. update, use prior update data
        estState.Pm = Pp
        estState.xm = xp

    # Get resulting estiamtes and variances
    posEst = estState.xm[:2]
    linVelEst = estState.xm[2:4]
    oriEst = estState.xm[4]
    driftEst = estState.xm[5]
    windEst = estState.xm[6]

    posVar = [estState.Pm[0,0], estState.Pm[1,1]]
    linVelVar = [estState.Pm[2,2], estState.Pm[3,3]]
    oriVar = estState.Pm[4,4]
    driftVar = estState.Pm[5,5]
    windVar = estState.Pm[6,6]

    return posEst,linVelEst,oriEst,windEst,driftEst,posVar,linVelVar,oriVar,windVar,driftVar,estState
