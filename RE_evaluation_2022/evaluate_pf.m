clear all
close all
clc

%%
delete(gcp('nocreate'));
euler = parcluster('local');
pool = parpool(euler, 48);

%%
addpath('PF_grading_2022');

nTrials = 20;

% Timeout per variation per student in seconds
timeoutS = 15 * 60;

% Error log
ferror = fopen('pf_error_log.txt', 'w+');

curDir = cd;
excelName = [curDir,'/pf_results.xls'];

submissions_dir = dir('Submissions');
submissions_dir = submissions_dir(3:end);
solutions_dir = dir('Solutions');
solutions_dir = solutions_dir(3:end);
D = [solutions_dir; submissions_dir];

EstimatorConstArray = {EstimatorConst(),EstimatorConstVar1(),EstimatorConstVar2()};
SimulationConstArray = {SimulationConst(),SimulationConstVar1(),SimulationConstVar2()};
n_variations = length(EstimatorConstArray);

data = cell(length(D)+1,1+4*(n_variations)+3);
header = ['Names', ...
          repmat({'Tracking Error','Nan Count','Avg. Ex. Time of thread pool', 'Avg. Ex. Time only the Estimator'}, 1, n_variations), ...
          'mean Tracking Error', 'mean Avg. Ex. Time of thread pool', 'mean Avg. Ex. Time only the Estimator'];
data(1, :) = header;

%% Run solutions and students code
for n = 1:length(D)
    
    estimator_path = fullfile(D(n).folder, D(n).name, 'P2_ParticleFilter');
    addpath(estimator_path);
    data(n + 1, 1) = {D(n).name};

    fprintf('Start evaluation #%d: %s\n', n, D(n).name)
    
    for var = 1:n_variations
        fprintf('Variation %d\n', var)

        SimConst = SimulationConstArray{var};
        EstConst = EstimatorConstArray{var};
                
        futures(1:nTrials) = parallel.FevalFuture;
        
        % Start stopwatch and evaluate runs asynchronously
        tic
        for i = 1:nTrials
            futures(i) = parfeval(pool, @run_estimator, 2, SimConst, EstConst, false, i);
        end

        % Wait for futures to complete and handle timeouts
        hasCompleted = wait(futures, 'finished', timeoutS);
        avg_time = toc/nTrials;
        
        trackingErrors = nan(nTrials, 1);
        tEstAvg = nan(nTrials, 1);

        if ~hasCompleted
            % Execution of runs has timed out without completing
            fprintf(ferror, ['Submission nr: %d\n' ...
                'Submission name: %s\n' ...
                'Error: timeout after %d seconds\n\n'], ...
                n, D(n).name, timeoutS);
            trackingErrors(1:nTrials) = nan;
        else
            % All trial runs have completed. Check for any potential errors
            for j = 1:nTrials
                if isempty(futures(j).Error)
                    % No error
                    [trackingErrors(j), tEstAvg(j)] = fetchOutputs(futures(j));
                else
                    % Error during run
                    fprintf(ferror, ['Submission nr: %d\n' ...
                        'Submission name: %s\n' ...
                        'Iteration number: %d\n' ...
                        'Error message: %s\n\n'], ...
                        n, D(n).name, j, ...
                        getReport(futures(j).Error, 'extended', 'hyperlinks', 'off'));
                    trackingErrors(j) = nan;
                end 
            end            
        end
        
        avg_trackingError = mean(trackingErrors, 'omitnan');
        avg_tEstAvg = mean(tEstAvg, 'omitnan');

        nancount = sum(isnan(trackingErrors));

        data(n+1,(var-1)*4+2:(var-1)*4+5) = {num2str(avg_trackingError), ...
            num2str(nancount),num2str(avg_time), num2str(avg_tEstAvg)};
        
    end

    %take the mean of n_variations and execution times
    data(n+1,1+4*n_variations+1) = {num2str(mean(str2double(convertCharsToStrings(data(n+1,2:4:10)))))};
    data(n+1,1+4*n_variations+2) = {num2str(mean(str2double(convertCharsToStrings(data(n+1,4:4:12)))))};
    data(n+1,1+4*n_variations+3) = {num2str(mean(str2double(convertCharsToStrings(data(n+1,5:4:13)))))};
    
    rmpath(estimator_path);
end


%%
rmpath('PF_grading_2022');
fclose(ferror);
save('pf_data');

writecell(data, 'pf_results.xls', 'UseExcel', false);

%%
%ws = load('data');
%data = ws.data;

%xlswrite(excelName, data, 1);

