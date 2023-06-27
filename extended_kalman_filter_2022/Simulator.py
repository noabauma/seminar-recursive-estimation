import numpy as np
import numpy.matlib                   #repmat
from scipy.integrate import solve_ivp #ode45 imitation

from SimulationConst import SimulationConst

""" helper functions """
def UniformCircle(R):
    r = np.random.rand()
    th = np.random.rand()
    x = np.sqrt(r)*np.cos(th)*R
    y = np.sqrt(r)*np.sin(th)*R
    return [x, y]

def Uniform(mx):
    return -mx + 2*mx*np.random.rand()

def UniformMinMax(mn,mx):
    return mn + (mx-mn)*np.random.rand()

def Normal(mx, vx):
    if vx < 0.0:
        raise Exception('Variance vx must be positive argument passed to "Normal"')
    
    return mx + np.sqrt(vx)*np.random.standard_normal()

def CalculateInputs(simConst: SimulationConst):
    # The length of the segment, in seconds.
    # Draw from uniform distribution.
    segTime = UniformMinMax(simConst.minSegTime,simConst.maxSegTime)

    # The applied inputs
    # Draw from uniform distribution. Note that the rudder angle is +-, 
    # while the thrust is only +.
    thrustForward = UniformMinMax(0,simConst.MaxThrust)
    rudderAngle = Uniform(simConst.MaxRudderCommand)

    return [thrustForward, rudderAngle, segTime]

def norm2(x):
    x = np.array(x)
    return x.dot(x)

def norm(x):
    x = np.array(x)
    return np.sqrt(x.dot(x))
# ----------------------------------------------------------------------- #



""" Simulation: Generate data """
# Simulation of the system to generate: inputs, measurements, true states.
# Inputs and measurements will be visible to the estimator.
#
# In this simulation, the inputs are piecewise constant functions of time, 
# that is for a sequence of random length, the input is constant.
def Simulator(simConst: SimulationConst):

    """ Initialize """
    # Number of samples of the simulation
    N = simConst.N

    # The starting position, velocity and orientation
    [posx_0, posy_0] =  UniformCircle(simConst.StartRadiusBound)
    vel_0 = np.zeros(2)
    angle_0 =  Uniform(simConst.RotationStartBound)

    # Starting wind direction
    wind_0 = Uniform(simConst.WindAngleStartBound)

    # The starting gyro drift
    drift_0 =  Uniform(simConst.GyroDriftStartBound)

    # The inputs
    inputThrust = np.zeros(N-1)
    inputAngle = np.zeros(N-1)

    # Time
    tm = np.round(np.linspace(0, (N-1)*simConst.sampleContinuous, N), 8) # TODO: python no round? Rounding due to floating point precision

    # The inputs for the first segment
    [thrustForward, rudderAngle, segTime] = CalculateInputs(simConst)

    # randomized sinusoid
    values = [-1.1, -1.0, -0.9, 0.9, 1.0, 1.1]
    sign_sin = values[np.random.randint(len(values))]
    offset = values[np.random.randint(len(values))]

    for n in range(N-1):
        # Store control input (known to estimator)
        inputThrust[n] = thrustForward
        inputAngle[n] = rudderAngle + 0.05*sign_sin*np.sin(tm[n]*0.01*2.0*np.pi) + 0.01*offset
        
        # Switch to new segment if we are at the end of the old one
        if tm[n] > segTime:
            [thrustForward, rudderAngle, segTime] = CalculateInputs(simConst)
            segTime = segTime + tm[n]

    input = np.array([inputThrust, inputAngle]) #this is for the return

    # Re-sample at 1kHz (well beyond the true dynamics)
    ssf = int(np.ceil(simConst.sampleContinuous/0.001))

    # Sample process noise
    # Constant 1000 because we sample at 1kHz -> adapt variance
    # n.shape = (ssf*N,4)
    n = np.hstack((np.sqrt(1000*simConst.DragNoise)*np.random.standard_normal(size=(ssf*N,1)),
                  np.sqrt(1000*simConst.RudderNoise)*np.random.standard_normal(size=(ssf*N,1)),
                  np.sqrt(1000*simConst.GyroDriftNoise)*np.random.standard_normal(size=(ssf*N,1)),
                  np.sqrt(1000*simConst.WindAngleNoise)*np.random.standard_normal(size=(ssf*N,1))))


    inputThrust = np.matlib.repmat(inputThrust, 1, ssf).T
    inputThrust = inputThrust.T.flatten()#.reshape(((N-1)*ssf,1))

    inputAngle = np.matlib.repmat(inputAngle, 1, ssf).T
    inputAngle = inputAngle.T.flatten()#.reshape(((N-1)*ssf,1))

    driftGyro = n[:, 2]

    driftWind = n[:, 3]

    # Simulate wind propagation
    wind_c = wind_0 + np.hstack((0 , np.cumsum(driftWind)/1000.0))
    wind   = wind_c[0:-1:ssf]


    # Simulate boat motion
    vdp1 = lambda t,x: [x[2],
                        x[3], 
                        np.cos(x[4])*(np.tanh(inputThrust[max(0,int(t/simConst.sampleContinuous*ssf))])-
                               simConst.dragCoefficientHydr*norm2([x[2], x[3]])*(1+n[max(0, int(t/simConst.sampleContinuous*ssf)),0]))-
                               simConst.dragCoefficientAir*(x[2]-simConst.windVel*np.cos(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))]))*
                               norm([x[2]-simConst.windVel*np.cos(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))]),
                               x[3]-simConst.windVel*np.sin(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))])]),
                        np.sin(x[4])*(np.tanh(inputThrust[max(0,int(t/simConst.sampleContinuous*ssf))])-
                               simConst.dragCoefficientHydr*norm2([x[2], x[3]])*(1+n[max(0, int(t/simConst.sampleContinuous*ssf)),0]))-
                               simConst.dragCoefficientAir*(x[3]-simConst.windVel*np.sin(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))]))*
                               norm([x[2]-simConst.windVel*np.cos(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))]),
                               x[3]-simConst.windVel*np.sin(wind_c[max(0,int(t/simConst.sampleContinuous*ssf))])]),
                        simConst.rudderCoefficient*inputAngle[max(0,int(t/simConst.sampleContinuous*ssf))]*(1+n[max(0,int(t/simConst.sampleContinuous*ssf)),1])]
    
    #sol = solve_ivp(vdp1, tm, [posx_0, posy_0, vel_0, angle_0])
    sol = solve_ivp(vdp1, [tm[0], tm[-1]], [posx_0, posy_0, vel_0[0], vel_0[1], angle_0], t_eval=tm)

    state = sol.y.transpose()

    # Simulate gyro drift
    drift_c = drift_0 + np.hstack((0, np.cumsum(driftGyro)/1000.0))
    drift = drift_c[0:-1:ssf]


    ## Measurements  
    # The distance sensors.  
    distASensor = np.zeros(N)
    distBSensor = np.zeros(N)
    distCSensor = np.zeros(N)

    # Sensor A,B & C measurement is acquired at randomized time intervals whose 
    # length is between the minimum and maximum length specified in the 
    # simulation constants. It is set to "np.inf" when no reading is made.

    # The Gyro sensor
    gyroSensor = np.zeros(N)

    # The Compass sensor
    compassSensor = np.zeros(N)

    # The next distance A,B & C reading: if the simulation time exceeds distCTime, a
    # measurement is acquired.
    distCTime =  UniformMinMax(simConst.sampleDistanceCmin,simConst.sampleDistanceCmax)

    for i in range(N):
        # The sensor reading.  The default is no measurement
        distASensor[i] = np.inf
        distBSensor[i] = np.inf
        distCSensor[i] = np.inf
        gyroSensor[i] = np.inf
        compassSensor[i] = np.inf

        
        
        # See if a distance A,B & C reading is made. 
        if tm[i] > distCTime:
            distCTime = tm[i] + UniformMinMax(simConst.sampleDistanceCmin,simConst.sampleDistanceCmax)

            distASensor[i] = norm(state[i,0:2] - simConst.pos_radioA) + Normal(0,simConst.DistNoiseA)
            distBSensor[i] = norm(state[i,0:2] - simConst.pos_radioB) + Normal(0,simConst.DistNoiseB)
            distCSensor[i] = norm(state[i,0:2] - simConst.pos_radioC) + Normal(0,simConst.DistNoiseC)
        
         
        
        gyroSensor[i] = state[i,4] + drift[i] + Normal(0,simConst.GyroNoise)
        compassSensor[i] = state[i,4] + Normal(0,simConst.CompassNoise)
    
    # No measurement for initial state (k=0)
    distASensor[0] = np.inf
    distBSensor[0] = np.inf
    distCSensor[0] = np.inf
    gyroSensor[0] = np.inf
    compassSensor[0] = np.inf
    
    sense = np.array([distASensor, distBSensor, distCSensor, gyroSensor, compassSensor])

    return tm, state, wind, drift, input.T, sense.T
