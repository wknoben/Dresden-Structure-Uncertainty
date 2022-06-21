% This file is part of the publication "Teaching hydrological modelling: 
% Illustrating model structure uncertainty with a ready-to-use teaching 
% module" (Knoben and Spieler, under review). It contains an example 
% of how the MARRMoT models m02 and m03 might be calibrated in the second
% part of the teaching module.
%
% NOTE: this file is derived from workflow_example_4.m that is part of the 
% MARRMoT v2.x distribution.
%
% Author:   Wouter J.M. Knoben
% Date:     22-06-2022
% Contact:  wouter.knoben@usask.ca
%
% Copyright (C) 2019, 2021, 2022 Wouter J.M. Knoben, Luca Trotter
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This example workflow includes 6 steps:
%
% 1. Data loading
% 2. Model choice
% 3. Model solver settings
% 4. Calibration settings
% 5. Model calibration within a loop over the two models and two catchments
% 6. Evaluation of calibration results

% NOTE: this example uses a custom function 'my_cmaes' to perform the
% optimisation, it is a wrapper around 'cmaes' to ensure inputs and outputs
% are consistent with other MATLAB's optimisation algorithms (e.g.
% 'fminsearch' or 'fminsearchbnd').
% While the wrapper is provided as part of MARRMoT, it requires the source 
% code to 'cmaes' to function, it is available at: 
% http://cma.gforge.inria.fr/cmaes.m
%
% The wrapper is necessary for the optimiser to function within the
% MARRMoT_model.calibrate method.
% Alternatively any model can be calibrated using any optimisation
% algorithm using the MARRMoT_model.calc_par_fitness method which returns
% the value of an objective function and can be used as input to an
% optimiser.

%% 1. Prepare data
% Load the data
load("Part 2 - catchment data.mat")

%% 2. Define model and catchments
models      = {'m_02_wetland_4p_1s','m_03_collie2_4p_1s'};
catchments  = {'c08109700','c12145500'};

%% 3. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver_opts.resnorm_tolerance = 0.1;                                  % Root-finding convergence tolerance; users have reported differences in simulation accuracy (KGE scores) during calibration between Matlab and Octave for a given tolerance. In certain cases, Octave seems to require tigther tolerances to obtain the same KGE scores as Matlab does.
input_solver_opts.rerun_maxiter   = 6;                                      % Maximum number of re-runs

%% 4. Define model-agnostic optimization settings
% Settings for 'my_cmaes'
% the opts struct is made up of two fields, the names are hardcoded, so
% they cannot be changed:
%    .sigma0:     initial value of sigma
%    .cmaes_opts: struct of options for cmaes, see cmaes documentation
%                 or type cmaes to see list of options and default values

% General options
optim_opts.TolFun       = 1e-6;                                             % stopping criterion on changes to fitness function
optim_opts.TolHistFun   = 1e-6;                                             % stopping criterion on changes to fitness function
optim_opts.EvalParallel = false;                                            % change to true to run in parallel on a pool of CPUs (e.g. on a cluster)

% Restarts?
% r = 2;                                                                    % uncomment to restart r-times until the condition
% optim_opts.Restarts  = r;
% optim_opts.RestartIf = 'fmin > -.8';                                      % OF is inverted, so this restarts unless the OF (KGE here) > 0.8

% Choose the objective function
of_name      = 'of_KGE';                                                    % This function is provided as part of MARRMoT. See ./MARRMoT/Functions/Objective functions
weights      = [1,1,1];                                                     % Weights for the three KGE components

% Time periods for calibration and evaluation. 
% Note: generally a 'warm-up period' is used to lessen the impact of the 
% initial conditions. Examples include running 1 year of data iteratively 
% until model stores reach an equilibrium, or choosing an arbitrary cut-off
% point before which the simulations are judged to be inaccurate and only 
% after which the objective function is calculated. For brevity, this step 
% is ignored in this example and the full calibration and evaluation 
% periods are used to calculate the objective function. Initial storages
% are estimated as 0 (see line 143 of this script).
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

% rewrite as CMAES inputs
n = length(time);
warmup = 0; 
cal_idx = (warmup+time_cal_start):time_cal_end;
eval_idx = (warmup+time_eval_start):time_eval_end;

% Remove all variables we no longer need
clear time*                                                                 % Remove date vector

%% 5. Calibrate the models
% MARRMoT model objects have a "calibrate" method that takes uses a chosen
% optimisation algorithm and objective function to optimise the parameter
% set. See MARRMoT_model class for details.

% Initiate the catchment loop
for c = 1:length(catchments)
    
    % Create a climatology data input structure. 
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_climatology.precip   = catchmentData.(catchments{c})(:,2);        % Daily data: P rate  [mm/d]
    input_climatology.temp     = catchmentData.(catchments{c})(:,3);        % Daily data: mean T  [degree C]
    input_climatology.pet      = catchmentData.(catchments{c})(:,4);        % Daily data: Ep rate [mm/d]
    input_climatology.delta_t  = 1;                                         % time step size of the inputs: 1 [d]
    
    % Extract observed streamflow
    Q_obs = catchmentData.(catchments{c})(:,5);                             % Daily data: Q rate [mm/d]
    
    % Initiate the model loop
    for im = 1:length(models)
    
        % Get the model settings and create the model object
        model     = models{im};                                             % Name of the model function (these can be found in Supporting Material 2)
        m         = feval(model);                                           % Model object
        parRanges = m.parRanges;                                            % Parameter ranges
        numParams = m.numParams;                                            % Number of parameters
        numStores = m.numStores;                                            % Number of stores
        input_s0  = zeros(numStores,1);                                     % Initial storages (see note in paragraph 5 on model warm-up)
        
        % Define the model-specific CMAES settings
        optim_opts.LBounds = parRanges(:,1);                                % lower bounds of parameters
        optim_opts.UBounds = parRanges(:,2);                                % upper bounds of parameters
        optim_opts.PopSize = 4 + floor(3*log(numParams));                   % population size (default)
        optim_opts.insigma = .3*(parRanges(:,2) - parRanges(:,1));          % starting sigma (this is default, could have left it blank)
        optim_opts.TolX    = 1e-6 * min(optim_opts.insigma);                % stopping criterion on changes to parameters 
        optim_opts.SaveFilename      = ['dresden_cmaesvars_',catchments{c},'_',models{im},'_2.mat']; % output file of cmaes variables
        optim_opts.LogFilenamePrefix = ['dresden_cmaeslogs_',catchments{c},'_',models{im},'_2.mat']; % prefix for cmaes log-files
        
        % Populate the model object
        m.input_climate = input_climatology;
        m.solver_opts   = input_solver_opts;
        m.S0            = input_s0;
        
        % Define an initial parameter set
        par_ini = 0.25*mean(parRanges,2);                                   % same as default value
        
        % Start calibration procedure
        disp('--- Calibration starting ---')
        
        % Check which objective function we start at
        Q_sim_ini = m.get_streamflow([],[],par_ini);
        of_start = feval(of_name,...                                        % Objective function name (here 'of_KGE')
                         Q_obs,...                                          % Observed flow during evaluation period
                         Q_sim_ini,...                                      % Simulated flow during evaluation period, using calibrated parameters            
                         cal_idx,...                                        % Indices of evaluation period
                         weights);                                          % KGE component weights
        disp(['Calibration starting at KGE = ',num2str(of_start)]) 
        
        % Call built-in calibration algorithm to optimize parameters
        [par_opt,...                                                        % optimal parameter set
         of_cal,...                                                         % value of objective function at par_opt
         stopflag,...                                                       % flag indicating reason the algorithm stopped
         output] = ...                                                      % other info about parametrisation
                    m.calibrate(...                                         % call the calibrate method of the model object
                                Q_obs,...                                   % observed streamflow
                                cal_idx,...                                 % timesteps to use for model calibration
                                'my_cmaes',...                              % function to use for optimisation (must have same structure as fminsearch)
                                par_ini,...                                 % initial parameter estimates
                                optim_opts,...                              % options to optim_fun
                                of_name,...                                 % name of objective function to use
                                1,...                                       % should the OF be inversed?
                                weights);                                   % additional arguments to of_name
                                   
        %% 6. Evaluate the calibrated parameters on unseen data
        % Run the model with calibrated parameters, get only the streamflow
        Q_sim = m.get_streamflow([],[],par_opt);
        
        % Compute evaluation performance
        of_eval = feval(of_name,...                                         % Objective function name (here 'of_KGE')
                        Q_obs,...                                           % Observed flow during evaluation period
                        Q_sim,...                                           % Simulated flow during evaluation period, using calibrated parameters            
                        eval_idx,...                                        % Indices of evaluation period
                        weights);                                           % KGE component weights
                           
        % Save the settings
        calibrationResults_v2.(catchments{c}).(models{im}).par      = par_opt;  % Calibrated parameter set
        calibrationResults_v2.(catchments{c}).(models{im}).KGE_cal  = of_cal;   % Calibration KGE
        calibrationResults_v2.(catchments{c}).(models{im}).KGE_eval = of_eval;  % Evaluation KGE
        
        % Print a summary
        Q_opt = m.get_streamflow([],[],par_opt);
        of_opt = feval(of_name,...                                          % Objective function name (here 'of_KGE')
                       Q_obs,...                                            % Observed flow during evaluation period
                       Q_opt,...                                            % Simulated flow during evaluation period, using calibrated parameters            
                       cal_idx,...                                          % Indices of evaluation period
                       weights);                                            % KGE component weights
        disp(['Calibration started at KGE = ',num2str(of_start)])
        disp(['Calibration ended at   KGE = ',num2str(of_opt)])
        disp(['Evaluation performance KGE = ',num2str(of_eval)])
        
    end
end

% Save results to disk
save calibrationResults_v2_clean.mat calibrationResults_v2