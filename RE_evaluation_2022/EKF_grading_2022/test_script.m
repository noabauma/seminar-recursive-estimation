close all;
clear;

%simConst = SimulationConst();
%[tm, state, wind, drift, input, sense] = Simulator(simConst);

N = 1000;

simConst = SimulationConst();
estConst = EstimatorConst();
doplot = false;

trackErrorNorm_arr = zeros(N,1);
angularErrorNorm_arr = zeros(N,1);
velocityErrorNorm_arr = zeros(N,1);
windErrorNorm_arr = zeros(N,1);
biasErrorNorm_arr = zeros(N,1);
time_arr = zeros(N,1);

for i = 1:N
    tic;
    [trackErrorNorm,angularErrorNorm,velocityErrorNorm,windErrorNorm,biasErrorNorm,initialEst,initialVar]=run_estimator(simConst,estConst,doplot,i);
    
    time_arr(i)              = toc;
    trackErrorNorm_arr(i)    = trackErrorNorm;
    angularErrorNorm_arr(i)  = angularErrorNorm;
    velocityErrorNorm_arr(i) = velocityErrorNorm;
    windErrorNorm_arr(i)     = windErrorNorm;
    biasErrorNorm_arr(i)     = biasErrorNorm;
end

%[tm, state, wind, drift, input, sense] = Simulator( simConst )

