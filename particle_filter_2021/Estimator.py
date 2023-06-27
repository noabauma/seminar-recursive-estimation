import numpy as np

from EstimatorConst import EstimatorConst

from PostParticles import PostParticles

"""
 The estimator function. The function will be called in two different
 modes: If km==0, the estimator is initialized. If km > 0, the
 estimator does an iteration for a np.single sample time interval unp.sing the 
 previous posterior particles passed to the estimator in
 prevPostParticles and the sensor measurement and control inputs.

 Inputs:
   prevPostParticles   previous posterior particles at time step k-1
                       The fields of the struct are [1xN_particles]-vector
                       (N_particles is the number of particles) 
                       corresponding to: 
                       .x_r: x-locations of the robot [m]
                       .y_r: y-locations of the robot [m]
                       .phi: headings of the robot [rad]
                       .kappa: wall offset [m]
                           
   sens                Sensor measurement z(k), scalar

   act                 Control inputs u(k-1), [1x2]-vector
                       act[0]: u_f, forward control input
                       act[1]: u_phi, angular control input

   estConst            estimator constants (as in EstimatorConst.m)

   km                  time index k, scalar
                       corresponds to continous time t = k*Ts
                       If km==0 initialization, otherwise estimator
                       iteration step.

 Outputs:
   postParticles       Posterior particles at time step k
                       The fields of the struct are [1xN_particles]-vector
                       (N_particles is the number of particles) 
                       corresponding to: 
                       .x_r: x-locations of the robot [m]
                       .y_r: y-locations of the robot [m]
                       .phi: headings of the robot [rad]
                       .kappa: wall offset [m]


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
"""
def Estimator(prevPostParticles: PostParticles, sens, act, estConst: EstimatorConst, km):

    # Set number of particles:
    N_particles = 1000 # obviously, you will need more particles than 10. But 10k is a bit slow...

    postParticles = PostParticles(3*N_particles)

    """ Mode 1: Initialization """
    if km == 0:
        # Do the initialization of your estimator here!
        n  = round(N_particles/2)
        ok = False
        while not ok:
            x = np.random.rand(2*N_particles)*2-1
            y = np.random.rand(2*N_particles)*2-1
            #accept = np.argwhere(np.np.sqrt(x**2 + y**2) <= 1.0)
            accept = np.argwhere(x**2 + y**2 <= 1.0)
            if len(accept) >= N_particles:
                ok = True
            
        
        
        postParticles.x_r = x[accept[:N_particles]]*estConst.d # N_particles matrix
        postParticles.y_r = y[accept[:N_particles]]*estConst.d # N_particles matrix
        postParticles.x_r = np.vstack((postParticles.x_r[:n]+estConst.pA[0],
                                       postParticles.x_r[n:]+estConst.pB[0])).reshape(-1)
        postParticles.y_r = np.vstack((postParticles.y_r[:n]+estConst.pA[1],
                                       postParticles.y_r[n:]+estConst.pB[1])).reshape(-1)
        
        postParticles.phi   = (np.random.rand(N_particles)*2-1)*estConst.phi_0 # N_particles matrix
        postParticles.kappa = (np.random.rand(N_particles)*2-1)*estConst.l # N_particles matrix
        
        # and leave the function
        return postParticles
     # end init

    """ Mode 2: Estimator iteration. """
    # If km > 0, we perform a regular update of the estimator.

    # Implement your estimator here!

    # Prior Update:
    v_f = (np.random.rand(N_particles)*2-1)*(estConst.sigma_f/2)
    v_phi = (np.random.rand(N_particles)*2-1)*(estConst.sigma_phi/2)
    x_prior = prevPostParticles.x_r + (act[0]+v_f) * np.cos(prevPostParticles.phi)
    y_prior = prevPostParticles.y_r + (act[0]+v_f) * np.sin(prevPostParticles.phi)
    phi_prior = prevPostParticles.phi + act[1] + v_phi
    k_prior = prevPostParticles.kappa

    # Posterior Update:

    #Find the distances by determining the intersection of the robot's direction with any boundary
    P=np.array([[0.5 ,0   ],
                [2.5 ,0   ],
                [2.5 ,1.5 ],
                [3   ,2   ],
                [2   ,3   ],
                [1.25,2.25],
                [1   ,2.5 ],
                [0   ,2.5 ], #x=kappa
                [0   ,0.5 ], #x=kappa
                [0.5 ,0.5 ]
                ])

    x_c = np.zeros(N_particles)
    y_c = np.zeros(N_particles)

    for n in range(N_particles):
        for i in range(10):
            A=np.array([[np.cos(phi_prior[n]), P[i,0]],
                        [np.sin(phi_prior[n]),P[i,1]]])
            if i == 9:
                A[:,1] = A[:,1] - P[0,:].T #closed room, need to reconnect P1
            else:
                A[:,1] = A[:,1] - P[i+1,:].T
            
            #handle the unknown kappa: harcoded, np.since the solution will not change
            if i==6:
                A[:,1] = [1-k_prior[n], 0]
            elif i==7:
                A[:,1] = [0, -2]
            elif i==8:
                A[:,1] = [0.5 - k_prior[n], 0]
            

            b = np.array([[P[i,0]-x_prior[n]],
                          [P[i,1]-y_prior[n]]])
            if i == 7 or i == 8:
                b[0] += k_prior[n]
            

            sol = np.linalg.solve(A,b)
            if sol[0] < 0.0: #intersection in the back of robot -> negative, no solution (parallel line) -> -Inf
                continue
            else:
                x_c[n] = np.round_(x_prior[n] + sol[0]*np.cos(phi_prior[n]),8)
                y_c[n] = np.round_(y_prior[n] + sol[0]*np.sin(phi_prior[n]),8)
                if i == 9:
                    if ((P[i,:] <= [x_c[n],y_c[n]]) & ([x_c[n],y_c[n]] <= P[0,:])).all():
                        #print("Heading towards line: " + i)
                        break
                    
                elif i < 3 or i == 8:
                    if ((P[i,:] <= [x_c[n],y_c[n]]) & ([x_c[n],y_c[n]] <= P[i+1,:])).all():
                        break
                    
                elif i == 3 or i == 5:
                    if ((P[i,0] >= x_c[n]) & (x_c[n] >= P[i+1,0]) & (P[i,1] <= y_c[n]) & (y_c[n] <= P[i+1,1])).all():
                        break
                    
                elif i == 4 or i == 6 or i == 7:
                    if ((P[i,:] >= [x_c[n],y_c[n]]) & ([x_c[n],y_c[n]] >= P[i+1,:])).all():
                        break
                    
                
            
    #         if i == 9:
    #             print("failed to assign line")
    #         
        
    
    #plt.scatter(x_c,y_c)

    #any deviation of estimated z[k] from sensed z[k] must be explained by w[k]
    w   = sens - np.sqrt((x_prior - x_c)**2 + (y_prior - y_c)**2)
    p_w = np.zeros(N_particles)

    #construct the function by region defintions:
    c1 = np.argwhere(abs(w) <= 2*estConst.epsilon)
    c2 = np.argwhere((2  *estConst.epsilon < abs(w)) & (abs(w) <= 2.5*estConst.epsilon))
    c3 = np.argwhere((2.5*estConst.epsilon < abs(w)) & (abs(w) <= 3  *estConst.epsilon))

    p_w[c1] = 2/(5*estConst.epsilon) - abs(w[c1])  /(5*estConst.epsilon**2)
    p_w[c2] =-4/(5*estConst.epsilon) + abs(w[c2])*2/(5*estConst.epsilon**2)
    p_w[c3] = 6/(5*estConst.epsilon) - abs(w[c3])*2/(5*estConst.epsilon**2)

    #threshold_p = 0.02 # if we np.argwhere the total probability below threshold, we resample anew
    if not any(p_w): #| sum(p_w)<threshold_p: # all particles have 0 (or very low) probability, i.e. we are lost
        #print("Lost Track, resample")
        ok=False
        while not ok: #check we sample inside the boundaries approximately
            postParticles.x_r = np.random.rand(3*N_particles)*3.2-0.2
            postParticles.y_r = np.random.rand(3*N_particles)*3
            resample1 = np.argwhere((postParticles.x_r < 0.5)  & (postParticles.y_r < 0.5)).reshape(-1)
            resample2 = np.argwhere((postParticles.x_r > 2.5)  & (postParticles.y_r < 1.5)).reshape(-1)
            resample3 = np.argwhere((postParticles.x_r < 1.25) & (postParticles.y_r > 2.5)).reshape(-1)
            resample4 = np.argwhere((postParticles.x_r + postParticles.y_r) > 5.0).reshape(-1)

            if len(resample1) + len(resample2) + len(resample3) + len(resample4) > N_particles:
                continue
            else:
                nosample = np.hstack((resample1,resample2,resample3,resample4))
                np.delete(postParticles.x_r, nosample)
                np.delete(postParticles.y_r, nosample)
                postParticles.x_r = postParticles.x_r[:N_particles]
                postParticles.y_r = postParticles.y_r[:N_particles]
                ok = True
            
        
        #postParticles.phi = np.random.rand(N_particles)*2*np.pi-np.pi
        postParticles.phi   = (np.random.rand(N_particles)*2-1)*np.pi
        postParticles.kappa = (np.random.rand(N_particles)*2-1)*estConst.l

        return postParticles #no further updates, just hope for better particle sampling this time
    

    alpha = 1.0/np.sum(p_w)
    beta = p_w*alpha

    betasum = beta
    for i in range(1,N_particles):
        betasum[i] += betasum[i-1]
    

    resampled_particles = np.zeros(N_particles, dtype=int)
    for i in range(N_particles):
        resampled_particles[i] = np.argwhere(betasum >= np.random.rand())[0]
    

    #print(size(unique(resampled_particles),1))

    #roughening
    d = 4.0
    K = 0.1
    E = np.array([[max(x_prior)  -min(x_prior)],
                  [max(y_prior)  -min(y_prior)],
                  [max(phi_prior)-min(phi_prior)],
                  [max(k_prior)  -min(k_prior)]]) #determine the maximal inter-sample variability
    
    sigmas=K*E*N_particles**(-1/d)

    postParticles.x_r   = x_prior[resampled_particles]   + np.random.standard_normal(N_particles)*sigmas[0]
    postParticles.y_r   = y_prior[resampled_particles]   + np.random.standard_normal(N_particles)*sigmas[1]
    postParticles.phi   = phi_prior[resampled_particles] + np.random.standard_normal(N_particles)*sigmas[2]
    postParticles.kappa = k_prior[resampled_particles]   + np.random.standard_normal(N_particles)*sigmas[3]

    while any(postParticles.phi > np.pi) or any(postParticles.phi < -np.pi):
        postParticles.phi[np.argwhere(postParticles.phi >  np.pi)]  = (np.random.rand()*2-1)*np.pi
        postParticles.phi[np.argwhere(postParticles.phi < -np.pi)] = (np.random.rand()*2-1)*np.pi
        #print("adjusted phi")
    

    return postParticles # end of estimator