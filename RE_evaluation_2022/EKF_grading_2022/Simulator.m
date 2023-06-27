function [tm, state, wind, drift, input, sense] = Simulator( simConst )
%% Simulation: Generate data
% Simulation of the system to generate: inputs, measurements, true states.
% Inputs and measurements will be visible to the estimator.
%
% In this simulation, the inputs are piecewise constant functions of time, 
% that is for a sequence of random length, the input is constant.

%% Initialize
% Number of samples of the simulation
N = simConst.N;

% The starting position, velocity and orientation
[posx_0, posy_0] =  UniformCircle(simConst.StartRadiusBound);
vel_0 = [0;0];
angle_0 =  Uniform(simConst.RotationStartBound);

% Starting wind direction
wind_0 = Uniform(simConst.WindAngleStartBound);

% The starting gyro drift
drift_0 =  Uniform(simConst.GyroDriftStartBound);

% The inputs
inputThrust = zeros(N-1,1);
inputAngle = zeros(N-1,1);

% Time
tm = round(linspace(0,(N-1)*simConst.sampleContinuous,N)',8); % Rounding due to floating point precision

% The inputs for the first segment
[thrustForward, rudderAngle, segTime] = CalculateInputs(simConst);

% randomized sinusoid
values = [-1.1 -1 -0.9 0.9 1 1.1];
sign_sin = values(randi(length(values)));
offset = values(randi(length(values)));

for n = 1:N-1
    % Store control input (known to estimator)
    inputThrust(n) = thrustForward;
    inputAngle(n) = rudderAngle + 0.05*sign_sin*sin(tm(n)*0.01*2*pi) + 0.01*offset;
    
    % Switch to new segment if we are at the end of the old one
    if (tm(n) > segTime)
        [thrustForward, rudderAngle, segTime] = CalculateInputs(simConst);
        segTime = segTime + tm(n);
    end    
end

input = [inputThrust inputAngle];

% Re-sample at 1kHz (well beyond the true dynamics)
ssf = ceil(simConst.sampleContinuous/0.001);

% Sample process noise
% Constant 1000 because we sample at 1kHz -> adapt variance
n = [sqrt(1000*simConst.DragNoise)*randn(ssf*N,1), ... 
     sqrt(1000*simConst.RudderNoise)*randn(ssf*N,1),...
     sqrt(1000*simConst.GyroDriftNoise)*randn(ssf*N,1),...
     sqrt(1000*simConst.WindAngleNoise)*randn(ssf*N,1)];

inputThrust = repmat(inputThrust,1,ssf)';
inputThrust = inputThrust(:);

inputAngle = repmat(inputAngle,1,ssf)';
inputAngle = inputAngle(:);

driftGyro = n(:,3);

driftWind = n(:, 4);

% Simulate wind propagation
wind_c = wind_0 + [0; cumsum(driftWind)/1000];
wind = wind_c(1:ssf:end-1);


% Simulate boat motion
[~,state] = ode45(...
    @(t,x)([x(3);...
            x(4);...
            cos(x(5))*(tanh(inputThrust(max(1,ceil(t/simConst.sampleContinuous*ssf))))-...
                simConst.dragCoefficientHydr*norm([x(3),x(4)])^2*(1+n(max(1,ceil(t/simConst.sampleContinuous*ssf)),1)))-...
                simConst.dragCoefficientAir*(x(3)-simConst.windVel*cos(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf)))))*...
                norm([x(3)-simConst.windVel*cos(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf)))),...
                x(4)-simConst.windVel*sin(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf))))]); ...
            sin(x(5))*(tanh(inputThrust(max(1,ceil(t/simConst.sampleContinuous*ssf))))-...
                simConst.dragCoefficientHydr*norm([x(3),x(4)])^2*(1+n(max(1,ceil(t/simConst.sampleContinuous*ssf)),1)))-...
                simConst.dragCoefficientAir*(x(4)-simConst.windVel*sin(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf)))))*...
                norm([x(3)-simConst.windVel*cos(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf)))),...
                x(4)-simConst.windVel*sin(wind_c(max(1,ceil(t/simConst.sampleContinuous*ssf))))]);...
            simConst.rudderCoefficient*inputAngle(max(1,ceil(t/simConst.sampleContinuous*ssf)))*(1+n(max(1,ceil(t/simConst.sampleContinuous*ssf)),2))]), ...
    tm, [posx_0;posy_0;vel_0;angle_0]);



% Simulate gyro drift
drift_c = drift_0 + [0; cumsum(driftGyro)/1000];
drift = drift_c(1:ssf:end-1);


%% Measurements  
% The distance sensors.  
distASensor = zeros(N,1);
distBSensor = zeros(N,1);
distCSensor = zeros(N,1);

% Sensor C measurement is acquired at randomized time intervals whose 
% length is between the minimum and maximum length specified in the 
% simulation constants. It is set to "inf" when no reading is made.

% The Gyro sensor
gyroSensor = zeros(N,1);

% The Compass sensor
compassSensor = zeros(N,1);

% The next distance A,B & C reading: if the simulation time exceeds distCTime, a
% measurement is acquired.
distCTime =  UniformMinMax(simConst.sampleDistanceCmin,simConst.sampleDistanceCmax);

for n = 1:N
    % The sensor reading.  The default is no measurement
    distASensor(n) = inf;
    distBSensor(n) = inf;
    distCSensor(n) = inf;
    gyroSensor(n) = inf;
    compassSensor(n) = inf;

    
       
    % See if a distance A,B & C reading is made. 
    if (tm(n) > distCTime)
        distCTime = tm(n) + UniformMinMax(simConst.sampleDistanceCmin,simConst.sampleDistanceCmax);
        
        distASensor(n) = norm(state(n,1:2) - simConst.pos_radioA) + Normal(0,simConst.DistNoiseA); %new measurements of A,B & C only taken when exceeding distCTime 
        distBSensor(n) = norm(state(n,1:2) - simConst.pos_radioB) + Normal(0,simConst.DistNoiseB); %new
        distCSensor(n) = norm(state(n,1:2) - simConst.pos_radioC) + Normal(0,simConst.DistNoiseC);
       
    end 
    
    gyroSensor(n) = state(n,5) + drift(n) + Normal(0,simConst.GyroNoise);
    compassSensor(n) = state(n,5) + Normal(0,simConst.CompassNoise);
end
% No measurement for initial state (k=0)
distASensor(1) = inf;
distBSensor(1) = inf;
distCSensor(1) = inf;
gyroSensor(1) = inf;
compassSensor(1) = inf;
 
sense = [distASensor, distBSensor, distCSensor, gyroSensor, compassSensor];

end


% ----------------------------------------------------------------------- %
% helper
function [x,y] = UniformCircle(R)
    r = rand;
    th = rand;
    x = sqrt(r)*cos(th)*R;
    y = sqrt(r)*sin(th)*R;
end

function val = Uniform(mx)
    val = -mx + 2*mx*rand;
end

function val = UniformMinMax(mn,mx)
    val = mn + (mx-mn)*rand;
end

function val = Normal(mx, vx)
    if(vx < 0)
        error('Variance vx must be positive argument passed to ''Normal''')
    end
    val = mx + sqrt(vx)*randn;
end

function [thrustForward, rudderAngle, segTime] = CalculateInputs(simConst)
    % The length of the segment, in seconds.
    % Draw from uniform distribution.
    segTime = UniformMinMax(simConst.minSegTime,simConst.maxSegTime);

    % The applied inputs
    % Draw from uniform distribution. Note that the rudder angle is +-, 
    % while the thrust is only +.
    thrustForward = UniformMinMax(0,simConst.MaxThrust);
    rudderAngle = Uniform(simConst.MaxRudderCommand);
end
