import numpy as np
from matplotlib import path
import copy

from EstimatorConst import EstimatorConst

from PostParticles import PostParticles

"""
 The estimator function. The function will be called in two different
 modes: If km==0, the estimator is initialized. If km > 0, the
 estimator does an iteration for a np.np.single sample time interval unp.np.sing the 
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

   km                  time index, scalar
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
                       .rho: wall offset for wall 1,2 [m]
                       .kappa: wall offset for wall 8,9 [m]


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
# ----------------------------------------------------------------------- #
# --------------------------helper functions----------------------------- #
# ----------------------------------------------------------------------- # 

# find closest point on contour
def findClosestPointOnContour(x,y,contour):
    dMin = np.inf
    for i in range(contour.shape[0]):
        if i < (contour.shape[0] - 1):
            dist, tmp1, tmp2 = lineToSegmentDistance(contour[i,0],contour[i+1,0],x,contour[i,1],contour[i+1,1],y)
        else:
            dist, tmp1, tmp2 = lineToSegmentDistance(contour[i,0],contour[0,0],x,contour[i,1],contour[0,1],y)

        if dist < dMin:
            contIdx = i
            dMin = dist

    # find coordinates onclosest contour
    if contIdx < contour.shape[0]-1:
        tmp, xC, yC = lineToSegmentDistance(contour[contIdx,0],contour[contIdx+1,0], x, contour[contIdx,1], contour[contIdx+1,1], y)
    else:
        tmp, xC, yC = lineToSegmentDistance(contour[contIdx,0],contour[0,0], x, contour[contIdx,1], contour[0,1], y)
        
    return xC, yC


# ------------------------------------------------------------------- #
def reInitEstimator(N_particles, estConst): 
    state = np.vstack((np.random.rand(N_particles)*2,
                       np.random.rand(N_particles)*2,
                       -np.pi + 2*np.pi*np.random.rand(N_particles),
                       UniformN(estConst.m, N_particles),
                       UniformN(estConst.l, N_particles) ))

    return state

# ------------------------------------------------------------------- #
def maxVariabilityAngles(angles):
    angles = myWrapTopi(angles)
    positive_angles = angles[angles >= 0]
    negative_angles = angles[angles < 0]
    
    diff_positive = max(positive_angles, default=0) - min(positive_angles, default=0)
    diff_negative = max(negative_angles, default=0) - min(negative_angles, default=0)

    tmp = np.vstack((np.hstack((negative_angles, positive_angles-np.pi)), np.hstack((np.zeros(negative_angles.shape[0]), np.ones(positive_angles.shape[0]) )) )).T
    
    mixed = tmp[tmp[:,0].argsort()] #same as sortrows
    
    temp = np.diff(mixed, axis=0)   #TODO: should be right? MATLAB gives a cool advice about the axis
    diff_mixed = min(abs(temp[0,temp[1,:] != 0]), default=0)

    return max([diff_positive, diff_negative, diff_mixed])     


# ------------------------------------------------------------------- #
def triTriangular(xBar, estConst):
    eps = estConst.epsilon
    if (0 <= abs(xBar)) and (abs(xBar) <= 2*eps):
        lHood = 2/(5*eps) - 1/(5*eps**2)*abs(xBar)            
    elif (2*eps < abs(xBar)) and (abs(xBar) <= 2.5*eps):
        lHood = 2/(5*eps**2)*(abs(xBar)-2*eps)
    elif (2.5*eps < abs(xBar)) and (abs(xBar) <= 3*eps):
        lHood = 1/(5*eps) - 2/(5*eps**2)*(abs(xBar)-5*eps/2) 
    else:
        lHood = 0
    
    return lHood

# ------------------------------------------------------------------- #
def propagatePosition(r,phi,u):
    return r + u*np.array([np.cos(phi),np.sin(phi)])

# ------------------------------------------------------------------- #
def getDist2Wall(xR,yR,phi,contour):
    Nwalls = contour.shape[0]
    wallsX = np.vstack((contour[:,0], np.hstack((contour[1:,0], contour[0,0])) )).T
    wallsY = np.vstack((contour[:,1], np.hstack((contour[1:,1], contour[0,1])) )).T
    distance2walls = np.inf*np.ones(Nwalls)

    for i in range(Nwalls):
        distance2walls[i] = computeRayLineIntersection(wallsX[i,0],wallsX[i,1],xR,wallsY[i,0],wallsY[i,1],yR,phi)

    return min(distance2walls)

# ------------------------------------------------------------------- #
def computeRayLineIntersection(x1, x2, xR, y1, y2, yR, phiR):
    w = [x1-xR, y1-yR]
    v = [np.cos(phiR), np.sin(phiR)]
    r = [x2-x1, y2-y1]
    # check if ray and line are parallel
    if abs(v[0] * r[1] - v[1] * r[0]) < 1e-10:
        return np.inf

    # t is the ray parameter; s is the line parameter
    t = (r[0] * w[1] - r[1] * w[0]) / (r[0] * v[1] - r[1] * v[0])
    s = (v[1] * w[0] - v[0] * w[1]) / (v[0] * r[1] - v[1] * r[0])

    # check if intersection is within line and on positive side of ray
    if t < 0 or s < 0 or s > 1:
        return np.inf
    else:
        return t
    

# ------------------------------------------------------------------- #
def lineToSegmentDistance(x1,x2,x3,y1,y2,y3):
    # (x1,y1) and (x2,y2) define the line, (x3,y3) is the point
    px = x2 - x1
    py = y2 - y1
    norm = px*px + py*py
    u =  ((x3 - x1) * px + (y3 - y1) * py) / norm
    if u > 1:
        u = 1
    elif u < 0:
        u = 0

    x = x1 + u * px
    y = y1 + u * py
    dx = x - x3
    dy = y - y3
    dist = np.sqrt(dx * dx + dy * dy)  # sqrt is not necessary for comparisons

    return dist, x, y

# ------------------------------------------------------------------- #

def myWrapTopi(x):
    x = np.array(x)
    q = (x < -np.pi) | (np.pi < x)
    x[q] = ((x[q] + np.pi) % (2*np.pi)) - np.pi

    return x


# ------------------------------------------------------------------- #
def UniformCircle(l, corners, N):
    # corners rows: [x0,y0]
    # number of rows: 2
    #samples = np.random.rand(3,N) # columns 1: corner, 2: radius, 3: angle
    samples = np.random.rand(3*N).reshape(N,3).T # for debugging
    r = np.sqrt(samples[1,:]) * l
    angle = samples[2,:] * 2 * np.pi
    corner1 = samples[0,:] > 0.5
    x = corner1 * corners[0,0] + (1.0 - corner1) * corners[1,0] + r * np.cos(angle)
    y = corner1 * corners[0,1] + (1.0 - corner1) * corners[1,1] + r * np.sin(angle)

    return x, y

# ----------------------------------------------------------------------- #
def Uniform(mx):
    return -mx + 2 * mx * np.random.rand()

# ----------------------------------------------------------------------- #
def UniformN(mx,N):
    return -mx + 2 * mx * np.random.rand(N)

def Estimator(prevPostParticles: PostParticles, sens, act, estConst: EstimatorConst, km):

    # Set number of particles:
    N_particles = 1000 # obviously, you will need more particles than 10.

    postParticles = PostParticles(3*N_particles)

    """ Mode 1: Initialization """
    if km == 0:
        # Do the initialization of your estimator here!
        
        # ################################################################### #
        # These particles are the posterior particles at discrete time k = 0
        # which will be fed into your estimator again at k = 1
        # Replace the following:
    
        corners = np.vstack((estConst.pA, estConst.pB))
        postParticles.x_r, postParticles.y_r = UniformCircle(estConst.d, corners, N_particles)
        postParticles.phi = UniformN(estConst.phi_0, N_particles)
        postParticles.rho = UniformN(estConst.m, N_particles)
        postParticles.kappa = UniformN(estConst.l, N_particles)

        # ################################################################### #
        # and leave the function
        return postParticles

    """ Mode 2: Estimator iteration. """
    # If km > 0, we perform a regular update of the estimator.

    # Implement your estimator here!
    # ####################################################################### #

    # Copy particles into np.single array:
    state = np.vstack((prevPostParticles.x_r, prevPostParticles.y_r, prevPostParticles.phi, prevPostParticles.rho, prevPostParticles.kappa ))

    """ Prior Update: """
    # Push particles through system dynamics:
    for p in range(N_particles):
        u0 = act[0] + Uniform(estConst.sigma_f/2)
        state[:2, p] = propagatePosition(state[:2,p].T, state[2,p], u0)
        state[2, p]  = myWrapTopi(state[2,p] + act[1] + Uniform(estConst.sigma_phi/2))

    """ Posterior Update: """
    # If there is a measurement, calculate the weights:
    betas = np.empty(N_particles)

    # Calculate expected measurement from prior particles
    for p in range(N_particles):
        contour_ = copy.deepcopy(estConst.contour)
        contour_[0,1] = state[3,p]
        contour_[1,1] = state[3,p]
        contour_[7,0] = state[4,p]
        contour_[8,0] = state[4,p]
        dist2Wall_p = getDist2Wall(state[0,p], state[1,p], state[2,p], contour_)
        betas[p] = triTriangular(sens - dist2Wall_p, estConst)

    # Redraw robot samples:
    cdf = np.cumsum(betas)
    stateTemp = np.zeros((5,N_particles))
    if cdf[-1] == 0:
        stateTemp = reInitEstimator(N_particles, estConst)
    else:
        for d in range(N_particles):
            ind = np.argwhere(cdf > np.random.rand()*cdf[-1])[0]
            stateTemp[:,d] = state[:,ind].ravel()


    state = stateTemp

    # Add roughening:
    Kr = 0.05
    D = 4
    # Calculate dmax, maximum inter-sample variability for each dimension:
    Emax = np.zeros(5)
    Emax[0] = np.max(state[0,:]) - np.min(state[0,:])
    Emax[1] = np.max(state[1,:]) - np.min(state[1,:])
    Emax[2] = maxVariabilityAngles(state[2,:])
    Emax[3] = np.max(state[3,:]) - np.min(state[3,:])
    Emax[4] = np.max(state[4,:]) - np.min(state[4,:])

    
    for dim in range(5):
        state[dim,:] += Kr * Emax[dim] * N_particles**(-1/D)# * np.random.normal(N_particles)


    # Ensure all particles are in valid range:
    contour_outside = copy.deepcopy(estConst.contour)
    contour_outside[0,1] = -estConst.m
    contour_outside[1,1] = -estConst.m
    contour_outside[7,0] = -estConst.l
    contour_outside[8,0] = -estConst.l

    polygon = path.Path(contour_outside)        #inpolygon
    inside = polygon.contains_points(state[:2,:].T)

    for idx in range(N_particles):
        if not inside[idx]:
            state[0,idx], state[1,idx] = findClosestPointOnContour(state[0,idx],state[1,idx],contour_outside)

    
    state[2,:] = myWrapTopi(state[2,:])

    postParticles.x_r = state[0,:]
    postParticles.y_r = state[1,:]
    postParticles.phi = state[2,:]
    postParticles.rho = state[3,:]
    postParticles.kappa = state[4,:]

        
    return postParticles # end of estimator