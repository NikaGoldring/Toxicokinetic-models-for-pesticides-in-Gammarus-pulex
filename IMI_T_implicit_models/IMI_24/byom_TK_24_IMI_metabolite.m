%% BYOM, byom_calanus_2016_onecomp.m 
%modified by Annika Mangold-D?ring 27.07.2020
%Intended use: Deriving TK parameter for minimilized TK of acute exposures
%(and TK experiment under different temperatures)
%
% * Author: Tjalling Jager 
% * Date: August 2018
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
global pri zvd      % global structures for optional priors and zero-variate data
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
%set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

DATA{1} = [ 1       17.57	17.57	17.57
            0.24	16.95	23.40	23.83
            1.01	69.93	74.25	64.62
            1.99	93.05	95.68	112.35
            2.94	68.56	92.55	68.71
            3.93	69.47	65.31	67.15
            5.08	55.05	65.18	57.53];

DATA{2} = [ 1	    17.57	17.57	17.57
            0.24	nan       nan      nan		
            1.01	nan       nan       nan		
            1.99	2.96	3.29	3.59
            2.94	3.04	6.06	5.03
            3.93	4.40	4.28	4.01
            5.08	3.01	4.59	4.21 ];


% W{1} = [10	9	10
%         9	10	10
%         10	9	10
%         9	7	8
%         9	7	9
%         10	8	9]; % each point is a pooled sample of X animals   % OR: 8 * ones(size(DATA{1})-1); 

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [17.57  % the scenarios (here nominal concentrations) 
         0
         0];  % initial values state 1 (internal concentrations)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.ke    = [0.2244    1 0.01 100 1];  % elimination rate constant, d-1
par.ku    = [3.723    1 0.01 100 1];  % uptake rate constant, L/kg/d
par.km    = [0.2    1 0.01 100 1];  % formation rate of the metabolite
par.kem   = [0.2    1 0.01 100 1];  % elimination rate of the metabolite

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['internal concentration (',char(181),'g/kg)'];
glo.ylab{2} = 'metabolite concentration';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (day)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'L']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

opt_optim.it = 0; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw  = 1; % plot in black and white
opt_plot.cn  = 1; % if set to 1, connect model line to points (only for bw=1)

glo.useode = 1; % AMD: 1 for using ODE solver, 0 for useing the analytical solution in simplefun.m

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the name of the parameter as it occurs in your parameter
% structure _par_ above. You do not need to run the entire script before
% you can make a profile. 
% 
% Options for the profiling can be set using opt_prof (see prelim_checks)

% Automatically calculate profiles for all parameters, and redo
% optimisation when a better value is found.
opt_prof.subopt = -1; % number of sub-optimisations to perform to increase robustness
% set to -1 to compare no sub-optims with 10 sub-optims.
opt_prof.detail = 1; % detailed (1) or a coarse (2) calculation

par_out = auto_profiles(par_out,opt_prof,opt_optim); 
% Here, the analysis showed that sub-optimisations are not needed

% print_par(par_out) % print best fit parameters in formatted manner

%% Likelihood region
% Another way to make intervals on model predictions is to use a sample of
% parameter sets taken from the joint likelihood-based conf. region. This
% is done by the function calc_likregion.m. It first does profiling of all
% fitted parameters to find the edges of the region. Then, Latin-Hypercube
% shooting, keeping only those parameter combinations that are not rejected
% at the 95% level in a lik.-rat. test. The inner rim will be used for CIs
% on forward predictions.
%
% Options for the likelihood region can be set using opt_likreg (see
% prelim_checks.m). 

opt_likreg.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_likreg.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_likreg.skipprof = 0; % skip profiling (and use profile from a saved mat file)
par_better = calc_likregion(par_out,500,opt_likreg); % second argument is target for number of samples in inner region (-1 to re-use saved sample from previous runs)

if isstruct(par_better) % if the profiling found a better optimum ...
    print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
    return % stop here, and don't go into plotting with CIs
end

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time

out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
