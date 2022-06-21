% This file is part of the publication "Teaching hydrological modelling: 
% Illustrating model structure uncertainty with a ready-to-use teaching 
% module" (Knoben and Spieler, under review). It contains an example 
% of how the MARRMoT models m02 and m03 might be calibrated in the second
% part of the teaching module.
%
% NOTE: this example uses a custom function 'fminsearchbnd', which is a
% basic optimization that lets the user specify constraints in the solution
% space. In this example these constraints are the model's parameter
% ranges. The file can be downloaded from Matlab's File Exchange:
% https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
%
% NOTE: this file is derived from workflow_example_4.m that is part of the 
% MARRMoT distribution.
%
% Author:   Wouter J.M. Knoben
% Date:     10-01-2021
% Contact:  wouter.knoben@usask.ca
%
% Copyright (C) 2021 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This example workflow includes 8 steps:
%
% 1. Check whether the optimization function is present
% 2. Data preparation
% 3. Model choice and setup
% 4. Model solver settings and time-stepping scheme
% 5. Calibration settings
% 6. Model calibration
% 7. Evaluation of calibration results

%% 1. Check function requirements
% Check whether fminsearchbnd is properly installed
if ~exist('fminsearchbnd','file')
    disp(['Optimizer ''fminsearchbnd'' not found. Visit the following ',...
          'link to download: ',...
          'https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon']);
    return
end

%% 2. Prepare data
% Load the data
load("Part 2 - catchment data.mat")

%% 3. Define model and catchments
models      = {'m_02_wetland_4p_1s','m_03_collie2_4p_1s'};
catchments  = {'c08109700','c12145500'};

%% 4. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver.name              = 'createOdeApprox_IE';                      % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;                                         % Maximum number of re-runs

%% 5. Define optimization settings
% Settings for 'fminsearchbnd'
optim_settings = optimset(...                                               % Using 'optimset' ensures that all non-specified fields are included in the resulting settings structure
                    'Display','iter',...                                    % Track progress (default = 'off')
                    'MaxIter',1000,...                                      % Stop after this many iterations (default = 200*numPar)
                    'MaxFunEvals',1000,...                                  % Stop after this many function evaluations (default = 200*numPar)
                    'TolFun',1e-4,...                                       % Stop if objective function change is below this tolerance (default = 1e-4)
                    'TolX',1e-4);                                           % Stop if changes in parameter values is below this tolerance (default = 1e-4)

% Choose the objective function
of_name      = 'of_KGE';                                                    % This function is provided as part of MARRMoT. See ./MARRMoT/Functions/Objective functions

% Time periods for calibration and evaluation. 
% Note: generally a 'warm-up period' is used to lessen the impact of the 
% initial conditions. Examples include running 1 year of data iteratively 
% until model stores reach an equilibrium, or choosing an arbitrary cut-off
% point before which the simulations are judged to be inaccurate and only 
% after which the objective function is calculated. For brevity, this step 
% is ignored in this example and the full calibration and evaluation 
% periods are used to calculate the objective function. Initial storages
% are estimated as 0 (see line 134 of this script).
%
% Calibration: 1989-01-01 to 1998-12-31
% Evaluation:  1999-01-01 to 2009-12-31
time            = datevec(catchmentData.(catchments{1})(:,1));              % Extract datenums as vectors
time_cal_start  = find(time(:,1) == 1989 & ...
                       time(:,2) == 01 & ...
                       time(:,3) == 01);                                    % Calibration start
time_cal_end    = find(time(:,1) == 1998 & ...
                       time(:,2) == 12 & ...
                       time(:,3) == 31);                                    % Calibration end
time_eval_start = find(time(:,1) == 1999 & ...
                       time(:,2) == 01 & ...
                       time(:,3) == 01);                                    % Evaluation start
time_eval_end   = find(time(:,1) == 2009 & ...
                       time(:,2) == 12 & ...
                       time(:,3) == 31);                                    % Evaluation end
clear time                                                                  % Remove date vector

%% 6. Calibrate the model
% Each MARRMoT model provides its outputs in a standardized way. We need
% the simulated flows for a given parameter set, so that we can compare
% simulated and observed flow, and determine the fitness of a particular
% parameter set. Simulated flow is stored in a structure of the form
% output.Q. The auxiliary function 'workflow_calibrationAssist' takes our
% specified model, inputs and a parameter set, runs the model and returns
% just simulated flow. The returned variable is used to optimize the
% parameter values. See lines 113-142.

% Initiate the catchment loop
for c = 1:length(catchments)
    
    % Create a climatology data input structure. 
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_climatology.precip   = catchmentData.(catchments{c})(:,2);        % Daily data: P rate  [mm/d]
    input_climatology.temp     = catchmentData.(catchments{c})(:,3);        % Daily data: mean T  [degree C]
    input_climatology.pet      = catchmentData.(catchments{c})(:,4);        % Daily data: Ep rate [mm/d]
    input_climatology.delta_t  = 1;                                         % time step size of the inputs: 1 [d]
    
    % Create temporary calibration time series
    input_climate_cal.precip  = input_climatology.precip(time_cal_start:time_cal_end);
    input_climate_cal.temp    = input_climatology.temp(time_cal_start:time_cal_end);
    input_climate_cal.pet     = input_climatology.pet(time_cal_start:time_cal_end);
    input_climate_cal.delta_t = input_climatology.delta_t;
    q_obs_cal                 = catchmentData.(catchments{c})(time_cal_start:time_cal_end,5);
    
    % Initiate the model loop
    for m = 1:length(models)
    
        % Model settings
        model    = models{m};                                               % Name of the model function (these can be found in Supporting Material 2)
        parRange = feval([model,'_parameter_ranges']);                      % Parameter ranges
        numPar   = size(parRange,1);                                        % Number of parameters
        numStore = str2double(model(end-1));                                % Number of stores
        input_s0 = zeros(numStore,1);                                       % Initial storages (see note in paragraph 5 on model warm-up)

        % Create a temporary function that calls the model, using the 
        % assisting function so that we get Qsim as an output, having only 
        % the parameter values as unknown input.
        q_sim_fun = @(par) workflow_calibrationAssist(...                   % Auxiliary function that returns Qsim
                                model,...                                   % Model name
                                input_climate_cal,...                       % Climate data during calibration period
                                input_s0,...                                % Initial storages
                                par,...                                     % We want to optimize these
                                input_solver);                              % Solver settings

        % Create the objective function we want to optimize. 
        % 'fminsearchbnd' is a minimizer and the objective function 
        % 'of_KGE' has range <-Inf,1] (<bad performance, good performance].
        % We thus want to optimize (minimize) -1*of_KGE which has range 
        % [-1,Inf> ([good,bad>).
        cal_fun = @(par) -1*(...                                                    
                            feval(of_name, ...                              % The objective function, here this is 'of_KGE'
                            q_obs_cal,...                                   % Observed flow during calibration period
                            q_sim_fun(par)));                               % Simulated flow for a given parameter set
                   
        % Create initial guesses for the optimizer
        par_ini = mean(parRange,2);
          
        % Call 'fminsearchbound' to optimize parameters
        disp('--- Calibration starting ---')
        [par_opt,of_cal] = fminsearchbnd(...                                % Returns optimized parameters and the objective function value
                            cal_fun,...                                     % Function that compares Qobs and Qsim for a given parameter set
                            par_ini,...                                     % Initial parameter guess
                            parRange(:,1),...                               % Lower parameter bounds
                            parRange(:,2),...                               % Upper parameter bounds
                            optim_settings);                                % Settings for fminsearchbnd
        
        %% 7. Evaluate the calibrated parameters on unseen data
        % Create temporary evaluation time series
        input_climate_eval.precip  = input_climatology.precip(time_eval_start:time_eval_end);
        input_climate_eval.temp    = input_climatology.temp(time_eval_start:time_eval_end);
        input_climate_eval.pet     = input_climatology.pet(time_eval_start:time_eval_end);
        input_climate_eval.delta_t = input_climatology.delta_t;
        q_obs_eval                 = catchmentData.(catchments{c})(time_eval_start:time_eval_end,5);
        
        % Run the model with calibration parameters
        model_out_eval = feval(model,...
                               input_climate_eval,...                       % Climate data during evaluation period
                               input_s0,...                                 % Initial storages
                               par_opt,...                                  % Calibrated parameters
                               input_solver);                               % Solver settings  
        
        % Compute evaluation performance
        of_eval = feval(of_name,...                                         % Objective function name (here 'of_KGE')
                        q_obs_eval,...                                      % Observed flow during evaluation period
                        model_out_eval.Q);                                  % Simulated flow during evaluation period, using calibrated parameters 
                           
        % Save the settings
        calibrationResults.(catchments{c}).(models{m}).par      = par_opt;  % Calibrated parameter set
        calibrationResults.(catchments{c}).(models{m}).KGE_cal  = -1*of_cal;% Calibration KGE, inverted back to <-Inf,1]
        calibrationResults.(catchments{c}).(models{m}).KGE_eval = of_eval;  % Evaluation KGE
        
    end
end

% Save results to disk
save calibrationResults.mat calibrationResults


%% Parameter uncertainty sampling
N = 10000; % Number of parameter samples to test

% Create empty storage structures
for c = 1:length(catchments)
    for m = 1:length(models)
        samplingResults.(catchments{c}).(models{m}).par = zeros(N,4);
        samplingResults.(catchments{c}).(models{m}).KGE_cal = zeros(N,1);
    end
end

% Initiate the catchment loop
for c = 1:length(catchments)
    
    % Create a climatology data input structure. 
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_climatology.precip   = catchmentData.(catchments{c})(:,2);        % Daily data: P rate  [mm/d]
    input_climatology.temp     = catchmentData.(catchments{c})(:,3);        % Daily data: mean T  [degree C]
    input_climatology.pet      = catchmentData.(catchments{c})(:,4);        % Daily data: Ep rate [mm/d]
    input_climatology.delta_t  = 1;                                         % time step size of the inputs: 1 [d]
    
    % Create temporary calibration time series
    input_climate_cal.precip  = input_climatology.precip(time_cal_start:time_cal_end);
    input_climate_cal.temp    = input_climatology.temp(time_cal_start:time_cal_end);
    input_climate_cal.pet     = input_climatology.pet(time_cal_start:time_cal_end);
    input_climate_cal.delta_t = input_climatology.delta_t;
    q_obs_cal                 = catchmentData.(catchments{c})(time_cal_start:time_cal_end,5);
    
    % Initiate the model loop
    for m = 1:length(models)
    
        % Model settings
        model    = models{m};                                               % Name of the model function (these can be found in Supporting Material 2)
        parRange = feval([model,'_parameter_ranges']);                      % Parameter ranges
        numPar   = size(parRange,1);                                        % Number of parameters
        numStore = str2double(model(end-1));                                % Number of stores
        input_s0 = zeros(numStore,1);                                       % Initial storages (see note in paragraph 5 on model warm-up)
        
        % Initiate sampling loop
        disp('--- Sampling starting ---')
        for n = 1:N
            
            % Get a random parameter sample
            par_trial = parRange(:,1) + (parRange(:,2)-parRange(:,1)) .* rand(numPar,1);
            
            % Run the trial parameter set to obtain simulated flows
            q_sim_trial = workflow_calibrationAssist(...                    % Auxiliary function that returns Qsim
                                model,...                                   % Model name
                                input_climate_cal,...                       % Climate data during calibration period
                                input_s0,...                                % Initial storages
                                par_trial,...                               % Randomly sampled parameter set
                                input_solver);                              % Solver settings
                            
            % Calculate the calibration KGE
            cal_kge = feval(of_name, ...                                    % The objective function, here this is 'of_KGE'
                            q_obs_cal,...                                   % Observed flow during calibration period
                            q_sim_trial);                                   % Simulated flow for a given parameter set
            
            % Save the parameter set and KGE score
            samplingResults.(catchments{c}).(models{m}).par(n,:) = par_trial;
            samplingResults.(catchments{c}).(models{m}).KGE_cal(n) = cal_kge;
        end
    end
end

% Save results to disk
save samplingResults.mat samplingResults 

%% Visualize
for c = 1:length(catchments)
    for m = 1:length(models)
        figure;
        for i=1:4 
            subplot(2,2,i); 
            scatter(samplingResults.(catchments{c}).(models{m}).par(:,i),samplingResults.(catchments{c}).(models{m}).KGE_cal,'.');
            ylim([1-sqrt(2),1])
            title([catchments{c}, ' ',models{m}])
        end
    end
end

