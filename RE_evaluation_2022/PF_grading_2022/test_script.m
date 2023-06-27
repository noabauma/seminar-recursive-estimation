clear all;
close all;
clc;

%rand('twister', 1);

simConst = SimulationConst();
estConst = EstimatorConst();

%simConst.N = 10;

%[km, state, input, sense] = Simulator(simConst);


n_runs = 1;

trackErrorNorm_arr = zeros(n_runs,1);
ExcTime_arr = zeros(n_runs,1);

for i = 1:n_runs
    tic;
    trackErrorNorm_arr(i) = run_estimator(simConst, estConst, true, 113);
    ExcTime_arr(i) = toc;
end