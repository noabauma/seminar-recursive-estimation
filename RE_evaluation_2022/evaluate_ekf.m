clear all
close all
clc

%%
delete(gcp('nocreate'));
euler = parcluster('local');
pool = parpool(euler, 48);

%%
addpath('EKF_grading_2022');

nTrials = 100;

% Timeout per variation per student in seconds
timeoutS = 15 * 60;

% Error log
ferror = fopen('ekf_error_log.txt', 'w+');

curDir = cd;
excelName = [curDir,'/ekf_results.xls'];

submissions_dir = dir('Submissions');
submissions_dir = submissions_dir(3:end);   %to ignore '.' & '..'
solutions_dir = dir('Solutions');
solutions_dir = solutions_dir(3:end);       %to ignore '.' & '..'
D = [solutions_dir; submissions_dir];

EstimatorConstArray = {EstimatorConst(),EstimatorConstVar1(),EstimatorConstVar2()};
SimulationConstArray = {SimulationConst(),SimulationConstVar1(),SimulationConstVar2()};
n_variations = length(EstimatorConstArray);

data = cell(length(D)+1,9*(n_variations+1)+2);
header = ['Names', ...
    repmat({'Tracking Error','Angular Error','Velocity Error','Wind Error','Bias Error','Initialization','Nan Count','Avg. Ex. Time of thread pool', 'Avg. Ex. Time only the Estimator'}, 1, n_variations), ...
    'mean & normalized Tracking Error','mean & normalized Angular Error','mean & normalized Velocity Error','mean & normalized Wind Error','mean & normalized Bias Error','mean Initialization','mean Nan Count','mean Avg. Ex. Time of thread pool', 'mean Avg. Ex. Time only the Estimator', 'overall error'];
data(1, :) = header;

%% Run solutions and students code

for n = 1:length(D)
    
    estimator_path = fullfile(D(n).folder, D(n).name, 'P1_EKF');
    addpath(estimator_path);
    data(n + 1, 1) = {D(n).name};

    fprintf('Start evaluation #%d: %s\n', n, D(n).name)
    
    for var = 1:n_variations
        fprintf('Variation %d\n', var)

        SimConst = SimulationConstArray{var};
        EstConst = EstimatorConstArray{var};
        
        trackingError = nan(nTrials,1);
        angularError = nan(nTrials,1);
        velocityError = nan(nTrials,1);
        windError = nan(nTrials,1);
        biasError = nan(nTrials,1);
        initialEst = nan(nTrials,6);
        initialVar = nan(nTrials,6);
        tEstAvg = nan(nTrials,1);
                
        futures(1:nTrials) = parallel.FevalFuture;
        
        % Start stopwatch and evaluate runs asynchronously
        tic
        for i = 1:nTrials
            futures(i) = parfeval(pool, @run_estimator, 8, SimConst, EstConst, false, i);
        end

        % Wait for futures to complete and handle timeouts
        hasCompleted = wait(futures, 'finished', timeoutS);
        avg_time = toc/nTrials;
        
        if ~hasCompleted
            % Execution of runs has timed out without completing
            fprintf(ferror, ['Submission nr: %d\n' ...
                'Submission name: %s\n' ...
                'Error: timeout after %d seconds\n\n'], ...
                n, D(n).name, timeoutS);
        else
            % All trial runs have completed. Check for any potential errors
            for j = 1:nTrials
                if isempty(futures(j).Error)
                    % No error
                    [trackingError(j),angularError(j),velocityError(j),windError(j), ...
                        biasError(j),initialEst(j,:),initialVar(j,:),tEstAvg(j)] = fetchOutputs(futures(j));
                else
                    % Error during run
                    fprintf(ferror, ['Submission nr: %d\n' ...
                        'Submission name: %s\n' ...
                        'Iteration number: %d\n' ...
                        'Error message: %s\n\n'], ...
                        n, D(n).name, j, ...
                        getReport(futures(j).Error, 'extended', 'hyperlinks', 'off'));
                end 
            end            
        end
        
        avg_trackingError = mean(trackingError, 'omitnan');
        avg_angularError = mean(angularError, 'omitnan');
        avg_velocityError = mean(velocityError, 'omitnan');
        avg_windError = mean(windError, 'omitnan');
        avg_biasError = mean(biasError, 'omitnan');
        
        if contains(D(n).folder, 'Solutions')
            initialEstTA = initialEst;
            initialVarTA = initialVar;
        end
        
        initPointsMean = length(find(sum(abs(initialEst - initialEstTA),1) < 1e-6));
        initPointsVar = length(find(sum(abs(initialVar - initialVarTA),1) < 1e-6));
        initPoints = initPointsMean+initPointsVar;

        nancount = sum(isnan(trackingError))+sum(isnan(angularError))+...
            sum(isnan(velocityError))+sum(isnan(biasError));
        
        avg_tEstAvg = mean(tEstAvg, 'omitnan');

        data(n+1,(var-1)*9+2:(var-1)*9+10) = {num2str(avg_trackingError),num2str(avg_angularError),...
            num2str(avg_velocityError),num2str(avg_windError),num2str(avg_biasError),...
            num2str(initPoints),num2str(nancount),num2str(avg_time), num2str(avg_tEstAvg)};
        
    end
    
    %take the mean of n_variations and execution times
    for i = 1:10
        data(n+1,9*n_variations+1+i) = {num2str(mean(str2double(convertCharsToStrings(data(n+1,1+i:9:9*n_variations+1)))))};
    end
    
    rmpath(estimator_path);
end

%normalize columns of the mean errors
for i = 29:33
   data(2:end,i) = cellstr(num2str(normalize(str2double(convertCharsToStrings(data(2:end,i))))));
end

%overall error is the sum of all normalized errors
for n = 1:length(D)
    data(n+1,end) = {num2str(sum(str2double(convertCharsToStrings(data(n+1,1+9*n_variations+1:1+9*n_variations+5)))))};
end

%%
rmpath('EKF_grading_2022');
fclose(ferror);
save('ekf_data');

writecell(data, 'ekf_results.xls', 'UseExcel', false);

%%
%ws = load('data');
%data = ws.data;

%xlswrite(excelName, data, 1);

