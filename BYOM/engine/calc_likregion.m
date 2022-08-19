function [par_better,CI_out] = calc_likregion(par,nr_lhs,opt_likreg)

% Usage: [par_better,CI_out] = calc_likregion(par,nr_lhs,opt_likreg)
%
% Calculates a likelihood-based joint confidence region. First, profile
% likelihoods for all fitted parameters are used to find the edges of the
% hyperbox that contains the true confidence region (of which the shape is
% then still unknown). Next, Latin Hypercube sampling is used to sample the
% hyperbox (shooting), and a likelihood-ratio test is used to decide which
% ones belong to the confidence region, and which are outside.
%
% This function also calculates the confidence intervals for the single
% parameters, just like <calc_proflik.m>. It has to do the profiling anyway
% to get the edges of the 'hyperbox', so you get the intervals on the
% single parameters for free. Second input is the number of samples that we
% aim for in the inner rim (within the chi2 criterion for df=1). This will
% be the set that is to be used for creating CIs on model predictions. Note
% that the we add a little bit to the chi2 criterion: since we use a
% discrete sample to approximate a 'hyper hull' in parameter space.
% 
% The likelihood profiling applies a variable stepsize algorithm. When the
% likelihood ratio changes very little (less than <Lcrit_min>), the
% stepsize is increased (up to a maximum, specified by <Fstep_max>). When
% the lik. ratio changes too much (more than <Lcrit_max>), the algorithm
% tries again with a smaller stepsize (also bound to a minimum:
% <Fstep_min>). Note that the stepsize is used as a fraction of the
% parameter value that is tried. To prevent very small stepsizes when the
% value goes towards zero (as can be the case for effect thresholds), I
% also apply an *absolute* minimum stepsize (<Fstep_abs>), which is
% specified as a fraction of the best parameter value (<Xhat>) (unless it
% is zero, then algoritm takes something small!).
%
% This function should run without the statistics toolbox of Matlab. But
% the LHS sampling will then be replaced by normal random sampling.
% Further, the critical value from the chi-square distribution will be read
% from a table, which goes to a maximum of 20 free parameters at the
% moment.
%
% For possible options to set in a structure as third argument
% (<opt_likreg>) see <prelim_checks.m>.
% 
% Author     : Tjalling Jager 
% Date       : March 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 h_txt

% read options from structure
prof_detail   = opt_likreg.detail; % detailed (1) or a coarse (2) calculation
if  ~ismember(prof_detail,[1 2])
    prof_detail = 2; % if a wrong value is entered, 'coarse' is default
end
n_sub        = opt_likreg.subopt; % number of sub-optimisations to perform to increase robustness
skip_prof    = opt_likreg.skipprof; % skip profiling step (1); use boundaries from saved set in MAT file from previous analysis
use_subplots = opt_likreg.subplt; % create single plot with subplots (1) or many individual plots (0)
chull        = opt_likreg.chull; % set to 1 to plot convex hull that approximates 95% edges of the likelihood region
axbnds       = opt_likreg.axbnds; % bind axes on the bounds of the hyperbox (1), accepted sample (2), or inner region (3)
burst        = opt_likreg.burst; % number of random samples from parameter space taken every iteration
brkprof      = opt_likreg.brkprof; % set to 1 to break the profiling when a better optimum is located
subann       = opt_likreg.subann; % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
lim_out      = opt_likreg.lim_out; % set to 1 to sample from a smaller part of space (enough for forward predictions)
Fsub_opt     = opt_likreg.subrng; % maximum factor on parameters (both higher and lower) for sub-optimisations

filenm       = glo.basenm;

par_better  = -1; % make sure the output is defined, even if we stop prematurely
names       = glo2.names;

sub_opt     = [Fsub_opt-1/Fsub_opt 1/Fsub_opt];% [2.67 0.33]; % settings for the suboptimisation: par * (rand * sub_opt(1)+ subopt(2))
sub_opt_log = [(log10(sum(sub_opt))-log10(sub_opt(2))) log10(sub_opt(2))]; % modify to get same range for log parameters!

% list of all critical values of the chi-square for 95% with df 1 to 20
% that way we can work without the statistics toolbox in most cases
chitable = [3.8415 5.9915 7.8147 9.4877 11.07 12.592 14.067 15.507 16.919 18.307 ...
    19.675 21.026 22.362 23.685 24.996 26.296 27.587 28.869 30.144 31.41];

diary ('results.out') % collect output in the diary "results.out"

disp(' ')
if exist('lhsdesign','file')~=2 % when lhsdesign does not exists as an m-file in the path
    warning('off','backtrace')
    Warning('You cannot use Latin-hypercube sampling; you need the statistics toolbox for that.')
    Warning('Instead, uniform random sampling will be used (more samples might be needed to obtain good coverage).')
    disp(' '), warning('on','backtrace')
end
if ~isfield(par,'tag_fitted') % apparently, parameters have not been fitted
    warning('off','backtrace')
    warning('Parameters have not been fitted, so results may not be very meaningful!')
    disp(' '), warning('on','backtrace')
end
    
names_tmp = names; % work with a copy of names, as the saved set may have different names!
% if the saved set has different names, this likely will produce an error
% elsewhere anyway.

pmat_tmp = packunpack(1,par,0); % use this to compare to the saved one, if re-using a saved set

if nr_lhs == -1 % do NOT make a new sample, but use the saved one!
    if exist([filenm,'_LR.mat'],'file') ~= 2
        error('There is no likelihood-region sample saved, so run calc_likregion again with a positive number of samples')
    end
    load([filenm,'_LR']) % load the random sample from the last likreg run
    % this loads the random sample in rnd, parameter matrix par and selection matrix par_sel
    % it also loads information from the profiling
    
    acc = rnd(rnd(:,end)<chitable(1),1:end-1); % keep the sets in inner rim (df=1) and remove last column
    rnd = rnd(:,1:end-1); % remove last column from rnd (likrat)
    
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
    if skip_prof == 0
        skip_prof = 1; % then also skip the profiling!
        disp('Skipping profiling as you requested to use the saved set.')
        disp(' ')
    end
    disp('Calculating sample from confidence region using previously determined bounds from MAT file. Profiles are reconstructed from the MAT file as well')
elseif skip_prof == 1 && exist([filenm,'_LR.mat'],'file') == 2
    disp('Skipping profiling; using bounds from saved set.')
    load([filenm,'_LR']) % load the profile information from the last calc_likregion run  
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
else
    if skip_prof == 1
        skip_prof = 0;
        disp('There is no saved set found, so calculating profiles anyway.')
        disp(' ')
    end
    disp('Calculating profiles and sample from confidence region ... please be patient.')
end
drawnow % plot the last plot in the plot buffer before starting the analysis

pmat = packunpack(1,par,0);  % transform structure into a regular matrix

if ~isequal(pmat,pmat_tmp)
    disp(' '), warning('off','backtrace')
    if isequal(pmat(:,1:2),pmat_tmp(:,1:2)) % ah, it's just the log settings or boundaries
        warning('The log settings and/or boundaries of parameters in the saved set differs from that in the workspace at the moment. This may not hinder the analysis.')
    else
        warning('The saved parameter matrix does not (exactly) equal the one entered when calling calc_likregion. The one from the saved set is used!')
    end
    warning('on','backtrace'), disp(' ')
    fprintf('Parameter values from saved set \n');
    fprintf('=================================================================================\n');
    nfields = length(names_tmp);
    if isfield(par,'tag_fitted')
        nfields = nfields - 1; % why is this needed?
    end
    for i = 1:nfields % display all parameters on screen
        if pmat(i,5) == 0
            fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g log-scale \n',names_tmp{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
        else
            fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g \n',names_tmp{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
        end
    end
    fprintf('=================================================================================\n');
    fprintf('  \n');
end

% put parameters that need to be on log scale on log10 scale (transfer needs it that way)
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
% also modify the min max ranges, but do that on a copy
pmat2 = pmat; % do this on a copy, as transfer wants the non-log version for the bounds
pmat2(pmat2(:,5)==0,[3 4]) = log10(pmat2(pmat2(:,5)==0,[3 4]));

parshat = pmat(:,1); % this is the best-fit parameter vector
par_sel = pmat(:,2); % this is the selection vector
ind_fit = find(par_sel == 1); % indext to fitted parameters
loglikmax = -1 * transfer(parshat(par_sel==1),pmat); % use transfer to obtain max log-likelihood

if length(ind_fit) == 1 % when we fit only 1 parameter ...
    use_subplots = 0; % turn of the subplots (no use for them)
end

% For finding the bounds of the joint hypercube region, use number of
% parameters as dfs. This is used at the bottom of this function to obtain
% a sample from the CI for prediction intervals. Making profiles here
% restricts the space that we need to sample.
if exist('chi2inv','file')==2 % when chi2inv exists exists as an m-file in the path ...
    chicrit  = chi2inv(0.95,sum(par_sel==1)); % to find the bounds of the hypercube, use the nr of estimated parameters as df
    chicrit2 = chi2inv(0.95,1); % to find 95% intervals for the single parameters and propagation region, use 1 df
else
    if sum(par_sel==1) > 20 % it is a bad idea to fit more than 20 parameters anyway
        error('For more than 20 fitted parameters, the critical value of the chi2 distribution was not included in calc_likregion.m.')
    end
    chicrit  = chitable(sum(par_sel==1)); % cut-off for the 95% parameter region
    chicrit2 = chitable(1); % cut-off for single par. CI and propagation region
end

if lim_out == 1
    chicrit = chicrit2 + 1; % just a bit more than the inner rim
end

% note that in pmat, log pars in first column are now on log scale!
pmat_orig = pmat;    % remember the original parameter values (we will be messing with it)
psel_orig = par_sel; % remember the original parameter selections (we will be messing with it)

%% Do the profiling

if use_subplots == 1
    [figh,ft] = make_fig(length(ind_fit),length(ind_fit)); % make figure of correct size
    hold on
end

sample_prof_acc = []; % initialise a matrix that will hold the points from the profile that will end up in the sample
% NOTE: the profiled, optimised, parameter sets will be added to the sample to increase robustness.

if skip_prof == 0

    boundscoll = zeros(length(ind_fit),2); % initialise the matrix to catch the bounds of the hyperbox
    fid = fopen('profiles_newopt.out','w'); % file to save all new optima as soon as the profile finds some
    % the 'w' option destroys the old contents of the file; use 'a' to append
    
    disp('  ')
    disp('Confidence interval from the profile (single parameter, 95% confidence)')
    disp('=================================================================================')
    
    flag_better  = 0; % keep track of whether a better value is found. 
    logliklowest = 0; % keep track of lowest minloglik found
    prof_coll = cell(length(ind_fit),1); % initialise cell array to catch profiles
    
    f = waitbar(0,['Calculating profile likelihoods for ',num2str(length(ind_fit)),' parameters. Please wait.'],'Name','calc_likregion.m');
    
    for i_p = 1:length(ind_fit)
        
        parnum  = ind_fit(i_p);   % find the index for fitted parameter i_p
        Xhat    = pmat(parnum,1); % the best estimate for the parameter (2nd column is select)
        
        % =====================================================================
        %               PERFORMANCE PARAMETERS
        % These values can be changed if the performance of the likelihood
        % procedure is not adequate.
        % =====================================================================
        
        switch prof_detail
            case 1 % detailed options
                Fstep_min = 0.001;    % min stepsize (as fraction of value that is tried)
                Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
                if Xhat == 0         % if parameter value is zero; than parameter is likely a NEC
                    Fstep_abs = 1e-4; % This is just a low value that should be ok in most situations
                else
                    Fstep_abs = 1e-3 * abs(Xhat); % smallest stepsize (in absolute sense)
                end
                % (this is to prevent problems when the parameter value that is tried goes
                % to zero; the stepsize is a fraction of the value tried so can also become
                % very small).
                
                % criteria for change in likelihood ratio
                Lcrit_max  = 1;   % maximum change in likelihood, above this value, stepsize is decreased
                Lcrit_min  = 0.2; % minimum change in likelihood, below this value, stepsize is increased
                Lcrit_stop = chicrit+5;  % stop when likelihood ratio reaches this value
                                
            case 2 % coarse options
                Fstep_min = 0.005;    % min stepsize (as fraction of value that is tried)
                Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
                if Xhat == 0         % if parameter value is zero; than parameter is likely a NEC
                    Fstep_abs = 1e-3; % This is just a low value that should be ok in most situations
                else
                    Fstep_abs = 1e-2 * abs(Xhat); % smallest stepsize (in absolute sense)
                end
                % (this is to prevent problems when the parameter value that is tried goes
                % to zero; the stepsize is a fraction of the value tried so can also become
                % very small).
                
                % criteria for change in likelihood ratio
                Lcrit_max  = 2; % maximum change in likelihood, above this value, stepsize is decreased
                Lcrit_min  = 0.7; % minimum change in likelihood, below this value, stepsize is increased
                Lcrit_stop = chicrit+2;  % stop when likelihood ratio reaches this value
                
        end
        
        % =====================================================================
        
        % depart from the original values (needed as we loop over all fitted parameters)
        par_sel = psel_orig;
        pmat    = pmat_orig; % note that log pars are on log scale in first column
        
        par_sel(parnum) = 0; % do not fit this parameter anymore, as we make profile
        pmat(:,2) = par_sel; % overwrite the selection part of pmat with the new par_sel
        
        p_index = (par_sel==1); % find the parameters that still need to be fitted
        Xhat = parshat(parnum); % best value of the profiled parameter
        
        % set options for the fitting (rougher than for finding the best fit; increases speed)
        oldopts  = optimset('fminsearch');
        options  = optimset(oldopts,'Display','off','FunValCheck','on'); % show no iterations
        options1 = optimset(options,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*(length(parshat)-1)); % only do rough optimisation
        options2 = optimset(options,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50*(length(parshat)-1)); % only do rough optimisation
        
        parmin  = pmat2(parnum,3); % use min and max of parameter as provided in script
        parmax  = pmat2(parnum,4); % this is from pmat2, so these may be on log-scale
        
        Xcoll     = []; % initialise the vector that collects the X-values
        loglikrat = []; % initialise the vector that collects the likrats
        
        for j = 1:2 % calculate from profile down and up from the ML estimate
            if j == 1
                % First calculating from max. lik. estimate downwards to zero ...
                proffact = -1; % going down ...
                
                if use_subplots == 0
                    figh2 = make_fig(1,1); % initialise a figure to keep track of the progress
                    h1 = gca; % remember the current axis number for plotting!
                    set(h1,'LineWidth',1,ft.name,ft.ticks); % adapt axis formatting
                    % if the parameter was on log-scale, mention it in the axis label
                    if pmat(parnum,5) == 0
                        xlabel(h1,['parameter ',names_tmp{parnum},' (log)'],ft.name,ft.label)
                    else
                        xlabel(h1,['parameter ',names_tmp{parnum}],ft.name,ft.label)
                    end
                    ylabel(h1,'2x loglik ratio','FontSize',12)
                    title(['Called from: ',filenm, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
                else 
                    figure(figh) % make sure that the multiplot is the current plot
                    g = subplot(length(ind_fit),length(ind_fit),(i_p-1)*length(ind_fit)+i_p); % subplot on diagonal
                    h1 = gca; % remember the current axis number for plotting!
                    set(h1,'LineWidth',1,ft.name,ft.ticks); % adapt axis formatting
                    p = get(g,'position');
                    p([3 4]) = p([3 4])*1.25; % add 10 percent to width and height
                    if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
                        p(2) = p(2)-0.04;
                    end
                    set(g, 'position', p);
                    g_rem(i_p) = g; % remember the handle to later correct axes if needed
                    
                    % determine which plots get a label on which axis
                    if i_p == 1
%                         if pmat(parnum,5) == 0
%                             ylabel(h1,[names_tmp{parnum},' (log)'],'FontSize',12)
%                         else
%                             ylabel(h1,[names_tmp{parnum}],'FontSize',12)
%                         end
                        ylabel(h1,'2x loglik ratio',ft.name,ft.label)
                        set(h1,'XTickLabel',[]); % this works with old versions as well
                    elseif i_p == length(ind_fit)
                        if pmat(parnum,5) == 0
                            xlabel(h1,[names_tmp{parnum},' (log)'],ft.name,ft.label)
                        else
                            xlabel(h1,[names_tmp{parnum}],ft.name,ft.label)
                        end
                        set(h1,'YTickLabel',[]); % remove tick labels on y-axis
                    else
                        set(h1,'YTickLabel',[]); % remove tick labels on y-axis
                        set(h1,'XTickLabel',[]); % remove tick labels on x-axis
                    end
                end
                
                hold on
                plot(h1,Xhat,0,'ko','MarkerFaceColor','r','MarkerSize',8) % plot the max lik value as a circle
                drawnow % update the plot
                             
            else
                % Now starting calculations from max. lik. estimate upwards ...
                proffact = 1; % going up ...
                % remember values from the first run (downwards) and clear 'em
                XcollD     = Xcoll;
                loglikratD = loglikrat;
                clear Xcoll loglikrat;
            end
            
            p_try     = pmat; % temporary parameter vector
            Xtry      = Xhat; % initialise the parameter to try on the ML estimate
            logliktmp = 0;
            i         = 1;
            
            % first values are zero
            Xcoll(1)     = Xtry; % all tried parameter values will be stored in Xcoll
            loglikrat(1) = 0;
            flag         = 0;
            
            % initialise Fstep as fraction of parameter value
            Fstep   = Fstep_min; % start from the minimum step size
            flagmin = 0; % test: flag is used to indicate that stepsize cannot or should not be decreased
            % this flag is put here to avoid trouble when the profile is followed
            % upwards. Each step, Fstep * Xtry will be higher than Fstep_abs, which
            % leads to a rejection of the value, a decrease in stepsize and a new
            % try.
            
            figure(f) % make sure that the waitbar is on top!
            
            while flag == 0 % as long as we have not reached a stopping criterion
                
                flag_try = 0;
                % remember the previous values
                Xtry_old      = Xtry;
                logliktmp_old = logliktmp;
                
                while flag_try == 0 && flag == 0 % as long as we haven't got a good new value
                    
                    if Xtry_old ==0 % if old tryvalue is zero, we cannot use the relative stepsize.
                        Xtry    = Xtry_old + proffact * Fstep_abs; % try a new value, based on absolute stepsize
                        Fstep   = Fstep_min; % and re-initialise the relative stepsize
                        flagmin = 1;
                    elseif Fstep * abs(Xtry_old) >= Fstep_abs % new stepsize is larger than abs min size
                        Xtry = Xtry_old + proffact * Fstep * abs(Xtry_old); % try a new value for the parameter
                        % use the given relative stepsize, or the minimum absolute stepsize
                        
                        % now, if we are at least 1.5 times above the absolute
                        % minimum and 1.5 times above the minimum relative stepsize
                        % ... then it makes sense to allow for a decrease in
                        % stepsize. Thus, flag is set to zero.
                        if Fstep * abs(Xtry_old) > 1.5 * Fstep_abs && Fstep > 1.5 * Fstep_min
                            flagmin = 0;
                        else
                            flagmin = 1;
                        end
                    elseif Fstep * abs(Xtry_old) < Fstep_abs
                        Xtry    = Xtry_old + proffact * Fstep_abs; % try a new value for the parameter
                        Fstep   = Fstep_abs / abs(Xtry_old);
                        flagmin = 1;
                    end
                    
                    % dont profile to negative values or very high ones!
                    % make sure that values are within the specified parameter bounds
                    % (bounds specified in calcstart.m)
                    Xtry = max(Xtry,parmin);
                    Xtry = min(Xtry,parmax);
                    
                    p_try(parnum,1) = Xtry ; % put new profile value into the total try vector
                    p_fit = p_try(p_index,1); % these are the ones that must be estimated now
                    
                    % Robustness of the optimisation routines is a real issue with
                    % profiling, unfortunately. Here, I implement two methods: use
                    % simulated annealling (probably less sensitive to local
                    % minima), or a method with simplex sub-optimisations.
                    
                    % minimisation is done with the reduced parameter set, other parameters in par are kept
                    % to the fixed value. This means that parshat will also contain the limited set only.
                    
                    if subann == 1
                        
                        opt_ann.Verbosity = 0; % controls output on screen (no output here)
                        opt_ann.StopTemp  = 1e-2; % crude temperature for stopping as Simplex will finish it off
                        [parshattemp,~] = anneal(@transfer,p_fit,opt_ann,p_try); % first rough estimation
                        
                    else
                        
                        % Now do fminsearch to get new best fit, with profile parameter fixed
                        % Do the optimisation at least twice: once rough, and once a little more
                        % precise. Hopefully, this will avoid local minima, but still
                        % be rapid enough!
                        
                        [parshatsub,logliksub] = fminsearch('transfer',p_fit,options1,p_try);
                        % Next, there is an option to perform additional
                        % sub-optimisations with perturbed starting values. This
                        % seriously increases robustness, but calculation time as well.
                        indmodlog = p_try(p_index,5)==0; % index to log parameters for fitted parameters
                        for i_sub = 1:n_sub
                            modpar = p_fit .* (rand(length(p_fit),1)*sub_opt(1)+sub_opt(2));
                            
                            if any(indmodlog==1) % for log-parameters, we need a different approach to get the same range
                                modpar(indmodlog) = p_fit(indmodlog) + rand(sum(indmodlog),1)*sub_opt_log(1)+sub_opt_log(2);
                            end
                            % make sure they are all within their bounds, otherwise fminsearch cannot start
                            modpar = max(modpar,pmat2(p_index,3));
                            modpar = min(modpar,pmat2(p_index,4));
                            
                            [parshatTMP2,logliknew] = fminsearch('transfer',modpar,options1,p_try);
                            % modify starting point by randomly making parameters higher and lower
                            logliksub(i_sub+1) = logliknew;
                            parshatsub(:,i_sub+1) = parshatTMP2;
                        end
                        [~,locsub] = min(logliksub); % find lowest loglik value from the sub-sampling
                        parshattemp = parshatsub(:,locsub); % and use the parameter values for that one
                        
                    end
                    
                    % use a Simplex to finish it
                    [parshattemp,logliknew] = fminsearch('transfer',parshattemp,options2,p_try);
                    
                    % and make the vector whole again
                    p_try(p_index,1) = parshattemp ; % replace values in total vector with values that were fitted
                    logliktmp        = -2 * (-1*logliknew-loglikmax); % likelihood ratio that follows a chi square
                    deltaL           = abs(logliktmp - logliktmp_old); % difference with previous likelihood
                    
                    if logliktmp <= chicrit % if the point is in the total 95% conf. region
                        sample_prof_acc  = [sample_prof_acc ;[p_try(ind_fit,1)' logliktmp]]; % add it as 'accepted set'
                    end
                    
                    % If the change in likelihood ratio is larger than we allow, or
                    % if we are at the border of the allowable parameter range, we
                    % want to decrease the stepsize. If the smalles stepsize is
                    % already used, we can stop the analysis.
                    skippit = 0;
                    
                    if deltaL > Lcrit_max
                        if flagmin == 0 % & logliktmp < 5,
                            % Only reject try and decrease stepsize when the flag says that it is possible!
                            % Test: and only when we are in an interesting range of likelihood values!
                            Fstep   = max(Fstep / 2 , Fstep_min); % try again with half the stepsize, unless less than min step
                            skippit = 1; % no storage for this one!
                        elseif any(Xtry == [parmin parmax])
                            flag = 1; % if already using smallest stepsize, stop!
                        end
                    elseif any(Xtry == [parmin parmax]) % also stop when max or min is reached and deltaL < Lcrit_max!
                        flag = 1; % stop!
                    end
                    
                    % If the change is small enough, or if we want to stop the analysis anyway, store the value.
                    if skippit == 0
                        flag_try = 1;
                        i = i+1;
                        loglikrat(i) = logliktmp; % this is 2x log-likelihood ratio
                        Xcoll(i)     = Xtry; % also remember the value of c0
                        plot(h1,Xcoll(i),loglikrat(i),'k.') % plot em on the fly!
                        drawnow
                        
                        if brkprof == 1 && flag_better == 1 && logliktmp > logliklowest
                            % We have found a better value than the best one
                            % already, but now it's getting worse again, so we can break this run.
                            par_better = packunpack(2,0,pmat_better); % return the better estimate in par_better
                            fprintf('  \n');
                            fprintf('Stopping the analysis: the profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,-1*loglikmax);
                            fprintf('Parameter values for new optimum \n');
                            fprintf('=================================================================================\n');
                            print_par(pmat_better); % write the better estimate to file in a formatted way
                            fprintf('=================================================================================\n');
                            fprintf('  \n');
                            diary off  % close results.out
                            return % so return to function calling us
                        end
                        
                        if logliktmp < logliklowest % if we found a lower likelihood, remember it!
                            logliklowest = logliktmp;
                            p_disp       = p_try;
                            loglik_disp  = logliknew;
                            
                            if logliktmp < -0.01
                                flag_better  = 1; % keep track of whether a better value is found.
                            end
                            
                            pmat_better  = [p_try(:,1) pmat_orig(:,2:5)]; % remember the better parameter set
                            % put parameters that were on log scale back on normal scale
                            pmat_better(pmat_better(:,5)==0,1) = 10.^(pmat_better(pmat_better(:,5)==0,1));
                            
                            % save it to the file profiles_newopt.out
                            ti = clock; % take current time
                            fprintf(fid,'%s (%02.0f:%02.0f) \n',date,ti(4),ti(5));
                            fprintf(fid,'The profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,-1*loglikmax);
                            fprintf(fid,'Parameter values for new optimum \n');
                            fprintf(fid,'=================================================================================\n');
                            print_par(pmat_better,fid); % write the better estimate to file in a formatted way
                            fprintf(fid,'=================================================================================\n');
                            fprintf(fid,'  \n');
                            
                        end
                        
                        if logliktmp > Lcrit_stop % stop if probability is high enough
                            flag = 1;
                        elseif deltaL < Lcrit_min % but change in likelihood is very small
                            Fstep = min(Fstep * 2 , Fstep_max); % so double the stepsize, unless greater than max step
                        end
                    end
                end
            end
            
            waitbar(((i_p-1)/length(ind_fit))+(0.5/length(ind_fit)),f) % make a nice waiting bar
            
        end
        
        waitbar(i_p/length(ind_fit),f) % make a nice waiting bar
        
        XcollU     = Xcoll; % the last Xcoll is for the run UPwards
        loglikratU = loglikrat;
        
        % ================================================
        % Now, derive the intercepts with the critical value for the joint conf. region
        % NOTE: this is improved by using the alternative method to
        % find zero crossings (see function calc_xing).
        
        prof_coll{i_p} = [XcollD' loglikratD';XcollU' loglikratU']; % construct entire profile
        prof_coll{i_p} = sortrows(prof_coll{i_p},1); % sort on basis of x-axis
        
        % find all crossings with the chi2 criterion
        XingS = calc_xing(prof_coll{i_p},chicrit2); % for single-par CI
        XingJ = calc_xing(prof_coll{i_p},chicrit);  % for total joint CR (unless lim_out=1)
        
        boundscoll(i_p,:) = [XingJ([1 end],1)]; % collect the extreme bounds for the confidence region
        % note that for the log-scale pars, boundscoll will stay on log scale!
        if use_subplots == 1 % set bound for the subplots at the confidence bounds
            xlim([boundscoll(i_p,1)-0.05*(boundscoll(i_p,2)-boundscoll(i_p,1)) boundscoll(i_p,2)+0.05*abs(boundscoll(i_p,2)-boundscoll(i_p,1))])
            ylim([min([loglikratD loglikratU]) chicrit+1])
        end
        
        % back-transform Xing2 if needed for CIs on screen
        if pmat(parnum,5) == 0
            XingS(:,1) = 10.^XingS(:,1);
        end 
        % column 2 is for warnings when hitting a boundary
        if XingS(1,2) == 1 % lowest crossing an upper boundary ...
            disp('Warning: taking the lowest parameter value as lower confidence limit')
        end
        if XingS(end,2) == 1 % highest crossing a lower boundary ...
            disp('Warning: taking the highest parameter value as upper confidence limit')
        end
        if size(XingS,1) == 2 % no problems, it is a single interval
            fprintf('%-10s interval: %#10.4g - %#1.4g \n',names_tmp{parnum},XingS(1,1),XingS(2,1))
        else
            disp('The confidence interval is a broken set (check likelihood profile to check these figures)')
            fprintf('%-10s interval: \n',names_tmp{parnum})
            for ix = 1:size(XingS,1)/2
                fprintf('   interval %1.1g: %#10.4g - %#1.4g \n',ix,XingS((ix-1)*2+1,1),XingS((ix-1)*2+2,1))
            end
        end
        
        CI_out(i_p,:) = [XingS(1,1) XingS(end,1)];
        
        % And make a nice plot
        plot(h1,prof_coll{i_p}(:,1),prof_coll{i_p}(:,2),'k-')
        
        % Plot the chi-square criterion for 95% and 1 degree of freedom
        plot(h1,[prof_coll{i_p}(1,1) prof_coll{i_p}(end,1)],[chicrit chicrit],'k:')
        plot(h1,[prof_coll{i_p}(1,1) prof_coll{i_p}(end,1)],[chicrit2 chicrit2],'k--')
        plot(h1,Xhat,0,'ko','MarkerFaceColor','r','MarkerSize',8) % plot the max lik value as a circle so its on top
        drawnow
        
        % % set the y-axis to allow to read the graph when lik-ratio becomes very large
        % plotmin = min(prof_coll{i_p}(:,2));
        % plotmax = max(10,min(30,max(prof_coll{i_p}(:,2)))); % force plotmax between 10 and 30
        
        if use_subplots == 0
            % ylim(h1,[plotmin plotmax])
            minX = prof_coll{i_p}(1,1);
            maxX = prof_coll{i_p}(end,1);
            xlim(h1,[minX-0.05*(maxX-minX) maxX+0.05*(maxX-minX)]) % limit x-axis
            if glo.saveplt > 0 % if we want to save the plot
                savenm = ['profile_reg_',filenm,'_par_',names_tmp{parnum}];%
                save_plot(figh2,savenm);
            end
        end
    end
    
    close(f) % close the waiting bar
    
    if lim_out == 0 % if we limit the space to search, chicrit is not the one for the joint CR
        disp('  ')
        disp('Edges of the joint 95% confidence region (using df=p)')
        disp('=================================================================================')
        
        for i_p = 1:length(ind_fit)
            parnum  = ind_fit(i_p);   % find the index for fitted parameter i_p
            if pmat(parnum,5) == 0 % parameter is fitted on log-scale
                fprintf('%-10s interval: %#10.4g - %#1.4g \n',names_tmp{parnum},10^boundscoll(i_p,1),10^boundscoll(i_p,2))
            else
                fprintf('%-10s interval: %#10.4g - %#1.4g \n',names_tmp{parnum},boundscoll(i_p,1),boundscoll(i_p,2))
            end
        end
        
        if logliklowest < 0 % if the profiler found a better optimum, display it!
            fprintf('\n  The profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,-1*loglikmax);
            % note that sampling will continue using the old loglikmax to
            % create the likelihood ratio ... perhaps better to use new best
            % ... although it is safest to start over from the new optimum
            fprintf('Parameter values for new optimum \n');
            fprintf('=================================================================================\n');
            print_par(pmat_better); % write the better estimate to screen in a formatted way
            fprintf('=================================================================================\n');
            fprintf('  \n');
        end
    end
    
    diary off    % close results.out
    
    % make the bounds a bit larger, just to be on the safe side
    boundscoll(:,1) = boundscoll(:,1)-0.05*(boundscoll(:,2)-boundscoll(:,1));
    boundscoll(:,2) = boundscoll(:,2)+0.05*(boundscoll(:,2)-boundscoll(:,1));
    
    fclose(fid); % close the file profiles_newopt.out for writing
    drawnow % plot the last plot in the plot buffer before doing the next analysis
    
    % Immediately save the profile information after profiling is
    % finished. This allows us to kill the next step if it takes too long,
    % and still have the profiles to start from.
    save([filenm,'_LR'],'par','par_sel','boundscoll','prof_coll','sample_prof_acc')
    % sample_prof_acc is added as these will also be included into the
    % random sample! Otherwise, these will be ignored if we redo the
    % analysis while skipping profiling.
    
elseif exist([filenm,'_LR.mat'],'file') == 2 && use_subplots == 1 % when we're not profiling, but there is a saved set ...
    % use that set to re-create the profiles, which makes the multiplot complete
    
    load([filenm,'_LR']) % load the random sample and profile info from the last calc_likregion run
    
    disp('  ')
    disp('Confidence interval from the profile (single parameter, 95% confidence)')
    disp('=================================================================================')
    
    for i_p = 1:length(ind_fit) % run through the fitted parameters
        figure(figh) % make sure that the multiplot is the current plot
        g = subplot(length(ind_fit),length(ind_fit),(i_p-1)*length(ind_fit)+i_p); % subplot on diagonal
        hold on
        h1 = gca; % remember the current axis number for plotting!
        set(h1,'LineWidth',1,'FontSize',10) % adapt axis formatting
        p = get(g,'position');
        p([3 4]) = p([3 4])*1.25; % add 10 percent to width and height
        if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
            p(2) = p(2)-0.04;
        end
        set(g, 'position', p);
        g_rem(i_p) = g; % remember the handle to later correct axes if needed
        
        % determine which plots get a label on which axis
        if i_p == 1
            if pmat(ind_fit(i_p),5) == 0
                ylabel(h1,[names_tmp{ind_fit(i_p)},' (log)'],'FontSize',12)
            else
                ylabel(h1,[names_tmp{ind_fit(i_p)}],'FontSize',12)
            end
            ylabel(h1,'2x loglik ratio','FontSize',12)
            set(h1,'XTickLabel',[]); % this works with old versions as well
        elseif i_p == length(ind_fit)
            if pmat(ind_fit(i_p),5) == 0
                xlabel(h1,[names_tmp{ind_fit(i_p)},' (log)'],'FontSize',12)
            else
                xlabel(h1,[names_tmp{ind_fit(i_p)}],'FontSize',12)
            end
            set(h1,'YTickLabel',[]); % remove tick labels on y-axis
        else
            set(h1,'YTickLabel',[]); % remove tick labels on y-axis
            set(h1,'XTickLabel',[]); % remove tick labels on x-axis
        end
        
        % And make a nice plot
        plot(h1,prof_coll{i_p}(:,1),prof_coll{i_p}(:,2),'k.-')
        
        % Plot the chi-square criterion for 95% and 1 degree of freedom
        if lim_out == 0 % if we limit the space to search, chicrit is not the one for the joint CR
            plot(h1,[min(prof_coll{i_p}(:,1)) max(prof_coll{i_p}(:,1))],[chicrit chicrit],'k:')
        end
        plot(h1,[min(prof_coll{i_p}(:,1)) max(prof_coll{i_p}(:,1))],[chicrit2 chicrit2],'k--')
        plot(h1,parshat(ind_fit(i_p)),0,'ko','MarkerFaceColor','r','MarkerSize',8) % plot the max lik value as a circle
        
        ylim([max(0,min(prof_coll{i_p}(:,2))) chicrit+1]) % limit y-axis to just above the cut-off criterion
        xlim([boundscoll(i_p,1) boundscoll(i_p,2)]) % limit x-axis to saved bounds
        drawnow
                
        % Display the CIs on screen as well
        parnum = ind_fit(i_p);   % find the index for fitted parameter i_p
        Xing   = calc_xing(prof_coll{i_p},chicrit2); % calculate crossing with a small function
        % back-transform Xing if needed
        if pmat(parnum,5) == 0
            Xing(:,1) = 10.^Xing(:,1);
        end 
        % column 2 is for warnings when hitting a boundary
        if Xing(1,2) == 1 % lowest crossing an upper boundary ...
            disp('Warning: taking the lowest parameter value as lower confidence limit')
        end
        if Xing(end,2) == 1 % highest crossing a lower boundary ...
            disp('Warning: taking the highest parameter value as upper confidence limit')
        end
        if size(Xing,1) == 2 % no problems, it is a single interval
            fprintf('%-10s interval: %#10.4g - %#1.4g \n',names_tmp{parnum},Xing(1,1),Xing(2,1))
        else
            disp('The confidence interval is a broken set (check likelihood profile to check these figures)')
            fprintf('%-10s interval: \n',names_tmp{parnum})
            for ix = 1:size(Xing,1)/2
                fprintf('   interval %1.1g: %#10.4g - %#1.4g \n',ix,Xing((ix-1)*2+1,1),Xing((ix-1)*2+2,1))
            end
        end
        
        CI_out(i_p,:) = [Xing(1,1) Xing(end,1)];
        
    end
    
end % this ends the skip_prof jump ...

%% Take an LHS (needs statistics toolbox) or normal random sample to find the joint conf. region
% pmat may contain parameters that are on log-scale, and boundscoll too.
% That is correct: transfer will put it back on normal scale.

if nr_lhs == -1 % do NOT make a new sample, but use the saved one!
    disp(' ')
    disp(['Number of samples accepted in the total conf. region: ',num2str(size(rnd,1)),', and in inner rim: ',num2str(size(acc,1))])
    disp(' ')
    
else
    % depart from the original values
    par_sel = psel_orig; % not used here, but is saved in the mat file
    pmat    = pmat_orig; % but parameter values in first column are on log scale

    nr_tot         = size(sample_prof_acc,1); % total number of parameter sets tried
    rnd = sample_prof_acc; % start from the profiled sets, and collect the accepted sets
    n_inner        = 0; % number of sets in inner rim (df=1)
    nr_it          = 0; % number of 'bursts'
    
    disp(' ')
    disp('Starting with obtaining a sample from the joint confidence region.')
    disp(['Bursts of ',num2str(burst),' samples, until at least ',num2str(nr_lhs),' samples are accepted in inner rim.'])
    
    f = waitbar(0,'Shooting for sample of joint confidence region. Please wait.','Name','calc_likregion.m');
    
    while n_inner < nr_lhs
        
        waitbar(n_inner/nr_lhs,f) % make a nice waiting bar

        if exist('lhsdesign','file')~=2 % when lhsdesign does not exist exists as an m-file in the path
            sample_lhs = rand(burst,length(ind_fit)); % uniform random sample between 0 and 1
        else
            sample_lhs = lhsdesign(burst,length(ind_fit)); % Latin-hypercube sample between 0 and 1
        end
        for i = 1:length(ind_fit) % go through the fitted parameters
            sample_lhs(:,i) = sample_lhs(:,i)*(boundscoll(i,2) - boundscoll(i,1))+boundscoll(i,1);
            % and change them to cover the bounds of the hypercube
        end
        loglik_lhs = zeros(size(sample_lhs,1),1); % pre-define for speed
        for i = 1:size(sample_lhs,1) % run through all samples
            loglik_lhs(i) = -1 * transfer(sample_lhs(i,:),pmat); % use transfer to obtain likelihood for each set
        end
        
        chi_lhs     = 2*(loglikmax - loglik_lhs); % calculate difference with the best fitting parameters that follows a chi-square
        ind_confreg = (chi_lhs<=chicrit); % this is the index to the sets that are within the region we like to keep
        rnd         = [rnd ;sample_lhs(ind_confreg,:) chi_lhs(ind_confreg)]; % only the accepted parameter values (within the joint search region)
        % also collect the minloglik ratio as last entry
        n_inner     = n_inner + sum(chi_lhs<=chicrit2); % how many are in the inner rim (df=1)
        
        % log-scale parameters are still on log-scale in this sample!
        nr_tot = nr_tot + size(sample_lhs,1);
        
    end
    
    close(f) % close the waiting bar
    
    % rnd is the accepted sample (within total 95% conf. region, df=p, or a limited region when lim_out = 1)
    rnd = sortrows(rnd,size(rnd,2)); % sort based on loglikrat (last column)
    
    % save the sample in a MAT-file with the name of the script, with
    % _LR at the end. Load can be used to retrieve it.
    % also, save par (total parameter matrix), par_sel (selection
    % matrix) and bounds to make plots without redoing the fit.
    save([filenm,'_LR'],'rnd','par','par_sel','boundscoll','prof_coll','sample_prof_acc')
        
    acc = rnd(rnd(:,end)<chicrit2,1:end-1); % put the sets in inner rim (df=1) into acc
    rnd = rnd(:,1:end-1); % remove likelihood-ratio column
    
    disp(['Size of total sample (profile points have been added): ',num2str(nr_tot)])
    disp(['of which accepted in the conf. region: ',num2str(size(rnd,1)),', and in inner rim: ',num2str(size(acc,1))])
    
end

%% Add plots with the sample in all binary combinations of parameters 

if axbnds == 1 % bind axis on the bounds of the hyperbox
    minrnd = boundscoll(:,1);
    maxrnd = boundscoll(:,2);
elseif axbnds == 2 % bind axis on the bounds of the 95% region
    minrnd = min(rnd,[],1);
    maxrnd = max(rnd,[],1);
elseif axbnds == 3 % bind axis on the bounds of the inner region
    minrnd = min(acc,[],1);
    maxrnd = max(acc,[],1);
end
% adding a bit extra on the bounds seems like a good idea in all cases
minrnd = minrnd - 0.05 * (maxrnd-minrnd);
maxrnd = maxrnd + 0.05 * (maxrnd-minrnd);
    
for parnr1 = 1:length(ind_fit)-1 % go through columns
    for parnr2 = parnr1+1:length(ind_fit) % go through rows
        
        if use_subplots == 0
            figh2 = make_fig(1,1); % initialise a new figure for every combination
            set(gca,'LineWidth',1,'FontSize',12,'FontName','Arial')
            title(['Called from: ',filenm, ' (',date,')'],'Interpreter','none','FontSize',10)
        else
            figure(figh) % make sure that the multiplot is the current plot
            g = subplot(length(ind_fit),length(ind_fit),(parnr2-1)*length(ind_fit)+parnr1);
            set(gca,'LineWidth',1,'FontSize',10,'FontName','Arial')
            p = get(g,'position');
            p([3 4]) = p([3 4])*1.25; % add to width and height
            if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
                p(2) = p(2)-0.04;
            end
            set(g, 'position', p);
        end
        ax = gca;
        hold on
        
        % Note that error ellipse is not relevant as it is not a random sample
        if lim_out == 0
            plot(rnd(:,parnr1),rnd(:,parnr2),'ko','MarkerFaceColor','c') % plot accepted set
        else
            plot(rnd(:,parnr1),rnd(:,parnr2),'co','MarkerFaceColor','w') % plot accepted set
        end
        plot(acc(:,parnr1),acc(:,parnr2),'ko','MarkerFaceColor','y') % plot sets within inner rim
        
        if chull == 1 && size(rnd,1)>3 % plot convex hull that approximates 95% edges of the likelihood region
            k = convhull(rnd(:,parnr1),rnd(:,parnr2));
            plot(rnd(k,parnr1),rnd(k,parnr2),'r-')
        end
        
        plot(parshat(ind_fit(parnr1)),parshat(ind_fit(parnr2)),'ko','MarkerFaceColor','r','MarkerSize',8) % plot best fit set
        axis([minrnd(parnr1) maxrnd(parnr1) minrnd(parnr2) maxrnd(parnr2)])
        
        if use_subplots == 0
            if pmat(ind_fit(parnr2),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                ylabel([names_tmp{ind_fit(parnr2)},' (log)'],'FontSize',12,'FontName','Arial')
            else
                ylabel([names_tmp{ind_fit(parnr2)}],'FontSize',12,'FontName','Arial')
            end
            if pmat(ind_fit(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                xlabel([names_tmp{ind_fit(parnr1)},' (log)'],'FontSize',12,'FontName','Arial')
            else
                xlabel([names_tmp{ind_fit(parnr1)}],'FontSize',12,'FontName','Arial')
            end
            
        else
            if parnr1 == 1 % we are in the first column
                if pmat(ind_fit(parnr2),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    ylabel([names_tmp{ind_fit(parnr2)},' (log)'],'FontSize',12,'FontName','Arial')
                else
                    ylabel([names_tmp{ind_fit(parnr2)}],'FontSize',12,'FontName','Arial')
                end
                if parnr2<length(ind_fit) % if we are not on the last row, remove axis numbers
                    set(ax,'XTickLabel',[]); % this works with old versions as well
                else
                    if pmat(ind_fit(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                        xlabel([names_tmp{ind_fit(parnr1)},' (log)'],'FontSize',12,'FontName','Arial')
                    else
                        xlabel([names_tmp{ind_fit(parnr1)}],'FontSize',12,'FontName','Arial')
                    end
                end
            elseif parnr2 == length(ind_fit) % we are in the last row
                if pmat(ind_fit(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    xlabel([names_tmp{ind_fit(parnr1)},' (log)'],'FontSize',12,'FontName','Arial')
                else
                    xlabel([names_tmp{ind_fit(parnr1)}],'FontSize',12,'FontName','Arial')
                end
                if parnr1>1 % if we are not on the first column, remove axis numbers
                    set(ax,'YTickLabel',[]); % this works with old versions as well
                end
            else
                set(ax,'YTickLabel',[]); % this works with old versions as well
                set(ax,'XTickLabel',[]); % this works with old versions as well
            end
        end
        
        if use_subplots == 0 && glo.saveplt > 0 % if we want to save the plot
            savenm = ['ellipse_reg_',filenm,'_',names_tmp{ind_fit(parnr2)},'x',names_tmp{ind_fit(parnr1)}];%
            save_plot(figh2,savenm);
        end
        
    end
end

if use_subplots == 1 % make sure the profiles have the same x-axis as the sample
    for i = 1:length(g_rem)
        xlim(g_rem(i),[minrnd(i) maxrnd(i)])
    end
end

if use_subplots == 1 % make an overall legend, top-right panel

    % create a legend in top-right sub-plot
    h_annot = subplot(length(ind_fit),length(ind_fit),length(ind_fit)); % make an extra subplot, and remember the handle
    p = get(h_annot,'position');
    p([3 4]) = p([3 4])*1.2; % add to width and height
    if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    h1 = gca; % remember the current axis
    set(h1,'LineWidth',1,ft.name,ft.ticks)
    dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
    hold on
    
    % plot fake data points, with same symbols as in regular plots
    if lim_out == 0
        plot(0,0,'ko','MarkerFaceColor','c') % plot total 95% conf. region
    else
        plot(0,0,'co','MarkerFaceColor','w') % plot limited sample
    end
    plot(0,0,'ko','MarkerFaceColor','y') % plot inner rim (df=1)
    plot(0,0,'ko','MarkerFaceColor','r','MarkerSize',8)
    if chull == 1 && size(rnd,1)>3 % plot convex hull that approximates 95% edges of the likelihood region
        plot([0 0],[0.1 0.1],'r-')
    end
    % And make a nice plot
    plot([0 0],[0.1 0.1],'k.-')
    
    % Plot the symbols for the profiles
    plot([0 0],[0.1 0.1],'k:')
    plot([0 0],[0.1 0.1],'k--')
    plot(0,0,'ko','MarkerFaceColor','r','MarkerSize',8) % plot the max lik value as a circle
    
    % define legend entries % removed:,['within limited set (',num2str(n_lim)' worst)']
    if lim_out == 0 % if we limit the space to search, chicrit is not the one for the joint CR
        L = {['within joint CI (df=',num2str(length(ind_fit)),')'], 'within inner rim (df=1)' ...
            'best fitting set','profile likelihood',['cutoff at df=',num2str(length(ind_fit))],'cutoff at df=1'};
    else
        L = {['within joint CI (limited)'], 'within inner rim (df=1)' ...
            'best fitting set','profile likelihood','within search region','cutoff at df=1'};
    end
    h_leg_tot = legend(L); % create a legend with all entries
    set(h_leg_tot,ft.name,ft.legend); % use standard font formats
    
    xlim([1 10]) % make sure the 'data' are off screen, so invisible
    dimleg = h_leg_tot.Position; % position of legend made
    % move it to the new subplot and set it to top-left of panel
    h_leg_tot.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];
    h1.Position = [dim(1) dim(2)+dim(4)-dimleg(4) 0.01 0.01];
    h1.Visible = 'off'; % hide axes
    
    if ~isfield(glo,'saveplt_notitle') || glo.saveplt_notitle ~= 1
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Likelihood-region method. Using file: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
    end
    % if use_subplots == 1 % || skip_prof == 1 % save the plot
        % if skip_prof == 1 && use_subplots == 0, figh = figh2, end % catch case when only making region, to save right plot
        if glo.saveplt > 0 % if we want to save the plot
            savenm = ['profile_reg_',filenm];%
            save_plot(figh,savenm);
        end
    % end
    
    
end

disp(' ')
disp(['Time required: ' secs2hms(toc)])
