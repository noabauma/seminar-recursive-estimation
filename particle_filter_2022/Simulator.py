import numpy as np

from SimulationConst import SimulationConst

# helper functions
two_pi = 2.0*np.pi
def wrapTo2Pi(x):
    tmp = x % two_pi
    positiveInput = x > 0.0
    tmp[(tmp == 0.0) & positiveInput] = two_pi #useless feature in my opinion
    return tmp

def wrapToPi(x):
    q = (x < -np.pi) | (np.pi < x)
    x[q] = wrapTo2Pi(x[q] + np.pi) - np.pi
    return x

def propagatePosition(r, phi, u):
    return r + u*np.array([np.cos(phi), np.sin(phi)])

# ----------------------------------------------------------------------- #
def UniformCircle(l, corners, N):
    # corners rows: [x0,y0]
    # number of rows: 2
    
    samples = np.random.rand(3,N) # rows 1: corner, 2: radius, 3: angle
    r = np.sqrt(samples[1,:])*l
    angle = samples[2,:]*2.0*np.pi

    corner1 = samples[0,:] > 0.5

    x = corner1*corners[0,0] + (1.0-corner1)*corners[1,0] + r*np.cos(angle)
    y = corner1*corners[0,1] + (1.0-corner1)*corners[1,1] + r*np.sin(angle)

    return x, y

# ------------------------------------------------------------------- #
def computeInput(r, phi, positive, simConst: SimulationConst):

    u_plus = UniformMinMax(0, simConst.alpha)
    u_minus = -u_plus

    if positive:
        r_plus = propagatePosition(r,phi,u_plus + simConst.sigma_f*0.5)
        d_plus = computeMinDistance(r_plus,simConst.contour)

        u = u_plus

        if d_plus < simConst.d_safety:
            r_minus = propagatePosition(r,phi,u_minus - simConst.sigma_f*0.5)
            d_minus = computeMinDistance(r_minus,simConst.contour)

            if d_minus > d_plus:
                u = u_minus
            elif d_minus == d_plus:
                u = 0.0
    else:
        r_minus = propagatePosition(r,phi,u_minus - simConst.sigma_f*0.5)
        d_minus = computeMinDistance(r_minus, simConst.contour)

        u = u_minus

        if d_minus < simConst.d_safety:
            r_plus = propagatePosition(r, phi, u_plus + simConst.sigma_f*0.5)
            d_plus = computeMinDistance(r_plus,simConst.contour)

            if d_plus > d_minus:
                u = u_plus
            elif d_minus == d_plus:
                u = 0.0

    return u

# ----------------------------------------------------------------------- #
def Uniform(mx):
    return -mx + 2.0*mx*np.random.rand()

# ----------------------------------------------------------------------- #
def UniformMinMax(mn,mx):
    return mn + (mx-mn) * np.random.rand()

# ----------------------------------------------------------------------- #
def getDistMeas(xR,yR,alphaR,contour,epsilon):
    # Compute distance measurement based on current state
    Nwalls = contour.shape[0]
    wallsX = np.vstack((contour[:,0], np.hstack((contour[1:,0], contour[0,0])) )).T
    wallsY = np.vstack((contour[:,1], np.hstack((contour[1:,1], contour[0,1])) )).T
    distance2walls = np.inf*np.ones(Nwalls)

    for i in range(Nwalls):
        distance2walls[i] = computeRayLineIntersection(wallsX[i,0],wallsX[i,1],xR,wallsY[i,0],wallsY[i,1],yR,alphaR)

    distance = min(distance2walls)

    # Add noise
    uSamp = np.random.rand()
    if uSamp <= 0.1:
        w_bar = inverseCumTrian(-3*epsilon,-2*epsilon,1/(5*epsilon),uSamp)
    elif (0.1 < uSamp) and (uSamp <= 0.9):
        uSamp -= 0.1
        w_bar = inverseCumTrian(-2*epsilon,2*epsilon,2/(5*epsilon),uSamp)

    elif 0.9 < uSamp:
        uSamp -= 0.9
        w_bar = inverseCumTrian(2*epsilon,3*epsilon,1/(5*epsilon),uSamp)

    return distance + w_bar

# ----------------------------------------------------------------------- #
def inverseCumTrian(a,b,c,y):
    x = 0.0                                             #Else is missing!? what is x else?
    if (0<=y) and (y <= c*(b-a)*0.25):
        x = a + np.sqrt(y*(b-a)/c)
    elif (c*(b-a)*0.25 <= y) and (y <= c*(b-a)*0.5):
        x = b - np.sqrt((b-a)**2*0.5 - y*(b-a)/c)       #<- bad code style: (b-a)/c*y
    return x

# ----------------------------------------------------------------------- #
def computeRayLineIntersection(x1,x2,xR,y1,y2,yR,alphaR):
    # helpers
    w = [x1-xR, y1-yR]
    v = [np.cos(alphaR), np.sin(alphaR)]
    r = [x2-x1, y2-y1]

    # check if ray and line are parallel
    if abs(v[0]*r[1]-v[1]*r[0]) < 1e-10:
        return np.inf

    # t is the ray parameter, s is the line parameter
    t = (r[0]*w[1]-r[1]*w[0])/(r[0]*v[1]-r[1]*v[0])
    s = (v[1]*w[0]-v[0]*w[1])/(v[0]*r[1]-v[1]*r[0])

    # check if intersection is within line and on positive side of ray
    if t < 0.0 or s < 0.0 or s > 1.0:
        return np.inf
    else:
        return t

# ----------------------------------------------------------------------- #
def computeMinDistance(r,contour):
    d = np.inf
    n = contour.shape[0]
    for i in range(n):

        if i < n-1:
            dist = lineToSegmentDistance(contour[i,0],contour[i+1,0],r[0],contour[i,1],contour[i+1,1],r[1])
        else:
            dist = lineToSegmentDistance(contour[i,0],contour[0,0],r[0],contour[i,1],contour[0,1],r[1])

        if dist < d:
            d = dist
    
    return d

# ----------------------------------------------------------------------- #
def lineToSegmentDistance(x1,x2,x3,y1,y2,y3):

    # (x1,y1) and (x2,y2) define the line, (x3,y3) is the point

    px = x2 - x1
    py = y2 - y1

    norm = px*px + py*py

    u =  ((x3 - x1) * px + (y3 - y1) * py) / norm

    if u > 1.0:
        u = 1.0
    elif u < 0.0:
        u = 0.0

    x = x1 + u * px
    y = y1 + u * py

    dx = x - x3
    dy = y - y3

    return np.sqrt(dx*dx + dy*dy) # sqrt is not necessary for comparisons

# ----------------------------------------------------------------------- #

def Simulator(simConst: SimulationConst):
    """ Simulation: Generate data"""
    # Simulation of the system to generate: inputs, measurements, true states.
    # Inputs and measurements will be visible to the estimator.
    
    simConst.alpha = 0.01    # Constant that defines the magnitude of the input
    simConst.d_safety = 0.1  # Safety distance from the wall

    """ Initialize """
    # Number of samples of the simulation
    N = simConst.N

    # The robot position
    r = np.zeros((N,2)) #[x_r,y_r]

    # The robot orientation
    phi = np.zeros(N)

    # Initialize uncertain points in contour
    rho = Uniform(simConst.m) * np.ones(N)
    simConst.contour[0,1] += rho[0]
    simConst.contour[1,1] += rho[0]

    kappa = Uniform(simConst.l) * np.ones(N)
    simConst.contour[7,0] += kappa[0]
    simConst.contour[8,0] += kappa[0]

    # The initial position and orientation
    corners  = np.vstack((simConst.pA, simConst.pB))
    x0, y0   = UniformCircle(simConst.d,corners,1)
    r[0,:]   = [x0,y0]
    phi[0] = Uniform(simConst.phi_0)

    # The input
    u = np.zeros((N-1,2)) #[u_f,u_phi]

    # The measurements
    distance = np.zeros(N-1)

    # Sign of last input
    positive = True

    for n in range(N-1):
        # Store control input (known to estimator)
        u[n,0] = computeInput(r[n,:], phi[n], positive, simConst)
        u[n,1] = 0.01
        
        # Update the last input sign
        positive = u[n,0] > 0.0
        
        # Add noise to the input
        u0 = u[n,0] + Uniform(simConst.sigma_f * 0.5)
    
        # Update position and orientation
        r[n+1,:] = propagatePosition(r[n,:], phi[n], u0)

        phi[n+1] = wrapToPi(np.array([phi[n] + u[n,1] + Uniform(simConst.sigma_phi * 0.5)]))

    
        # Store distance measurement    
        distance[n] = getDistMeas(r[n+1,0], r[n+1,1], phi[n+1], simConst.contour, simConst.epsilon)

    km = np.arange(simConst.N+1)#.reshape((-1,1))
    
    phi = phi.reshape((-1,1))
    rho = rho.reshape((-1,1))
    kappa = kappa.reshape((-1,1))
    state = np.hstack((r, phi, rho, kappa))
    input = u
    sense = np.hstack((np.inf, distance)).reshape((-1,1))

    return km, state, input, sense
# ----------------------------------------------------------------------- #