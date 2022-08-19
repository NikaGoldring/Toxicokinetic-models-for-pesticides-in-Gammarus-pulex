function [Xing,par_better,prof_coll,loglik_disp] = calc_proflik(par,parname,opt_prof)

% Usage: [Xing,par_better,prof_coll,loglik_disp] = calc_proflik(par,parname,opt_prof)
%
% Calculates a profile likelihood for parameter with the name <parname>
% from your parameter structure. <par> is the total parameter structure,
% resulting from the optimisation. Various options can be set in <opt_prof>
% (see <prelim_checks.m>), a.o. to specify the level of detail of the
% profiling, detailed (1) or coarse (2, default). Optional outputs: edges
% of the 95% CI in <Xing>, better parameter estimates in <par_better>
% (might be only very slightly better), entire profile in <prof_coll> (log
% pars on log scale), and better log-likelihood if a better value was
% found.
%
% This function is slightly modified from the function in DEBtoxM.
%
% This function applies a variable stepsize algorithm. When the likelihood
% ratio changes very little (less than <Lcrit_min>), the stepsize is
% increased (up to a maximum, specified by <Fstep_max>). When the lik.
% ratio changes too much (more than <Lcrit_max>), the algorithm tries again
% with a smaller stepsize (also bound to a minimum: <Fstep_min>). Note that
% the stepsize is used as a fraction of the parameter value that is tried.
% To prevent very small stepsizes when the value goes towards zero (as can
% be the case for effect thresholds), I also apply an *absolute* minimum
% stepsize (<Fstep_abs>), which is specified as a fraction of the best
% parameter value (<Xhat>) (unless it is zero, then algoritm takes
% something small!).
%
% Author     : Tjalling Jager
% Date       : April 2019
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2

names = glo2.names;

% predefine outputs, so they are known if function returns prematurely
par_better = -1;
Xing       = [];
prof_coll  = -1;

% read options from structure
prof_detail = opt_prof.detail; % detailed (1) or a coarse (2) calculation
if  ~ismember(prof_detail,[1 2])
    prof_detail = 2; % if a wrong value is entered, 'coarse' is default
end
n_sub         = opt_prof.subopt;  % number of sub-optimisations to perform to increase robustness
brkprof       = opt_prof.brkprof; % set to 1 to break the profiling when a better optimum is located
verbose       = opt_prof.verbose; % set to 0 to suppress all output from calc_proflik to screen
subann        = opt_prof.subann;  % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
Fsub_opt      = opt_prof.subrng; % maximum factor on parameters (both higher and lower) for sub-optimisations
scriptnm      = glo.basenm;

sub_opt     = [Fsub_opt-1/Fsub_opt 1/Fsub_opt];% [2.67 0.33]; % settings for the suboptimisation: par * (rand * sub_opt(1)+ subopt(2))
sub_opt_log = [(log10(sum(sub_opt))-log10(sub_opt(2))) log10(sub_opt(2))]; % modify to get same range for log parameters!

parnum = find(strcmp(names,parname)); % lookup the name in the structure and return the position
if isempty(parnum)
    error(['The name of the parameter to profile (',parname,') does not exist in your structure.'])
end
if ~isfield(par,'tag_fitted') % apparently, parameters have not been fitted
    warning('off','backtrace')
    warning('Parameters have not been fitted, so results may not be very meaningful!')
    disp(' '), warning('on','backtrace')
end

fid = fopen('profiles_newopt.out','w'); % file to save all new optima as soon as the profile finds some
% the 'w' option destroys the old contents of the file; use 'a' to append 

pmat    = packunpack(1,par,0);  % transform structure into a regular matrix
if size(pmat,2) < 4 % no full parameter definition with bounds
    error('Specify the parameters in your script file with upper and lower bounds.')
end
if size(pmat,2) == 4 % when not fitted first, the 5th column is not added if needed
    pmat(:,5) = 1; % by default, assume normal scale
end
% put parameters that need to be on log scale on log10 scale
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
% also modify the min max ranges!
pmat2 = pmat; % but do this on a copy, as transfer wants the non-log version of the ranges
pmat2(pmat2(:,5)==0,[3 4]) = log10(pmat2(pmat2(:,5)==0,[3 4]));

parshat = pmat(:,1);            % this is the parameter vector
par_sel = pmat(:,2);            % this is the selection vector
Xhat    = pmat(parnum,1);       % the best estimate for the parameter (2nd column is select)

%% Set performance parameters
% These values can be changed if the performance of the likelihood
% procedure is not adequate.

chicrit = 3.8415; % = chi2inv(0.95,1); % the chi square criterion for 1 df and 95%
% This is coded in hard form so that this function can operate without the
% statistics toolbox.

switch prof_detail
    case 1 % detailed options
        Fstep_min = 0.001;    % min stepsize (as fraction of value that is tried)
        Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
        if Xhat == 0          % if parameter value is zero; than parameter is likely a NEC
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
        Lcrit_stop = 15;  % stop when likelihood ratio reaches this value
        
    case 2 % coarse options
        Fstep_min = 0.005;    % min stepsize (as fraction of value that is tried)
        Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
        if Xhat == 0          % if parameter value is zero; than parameter is likely a NEC
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
        Lcrit_stop = 10;  % stop when likelihood ratio reaches this value
        
end

% set options for the fitting (rougher than for finding the best fit; increases speed)
oldopts  = optimset('fminsearch');
options  = optimset(oldopts,'Display','off','FunValCheck','on'); % show no iterations
options1 = optimset(options,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*(length(parshat)-1)); % only do rough optimisation
options2 = optimset(options,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50*(length(parshat)-1)); % only do rough optimisation

%% Start with the profiling ...

par_sel(parnum) = 0; % do not fit this parameter anymore, as we make profile
pmat_orig = pmat;    % remember the original parameter values (we are messing with a copy)
pmat(:,2) = par_sel; % overwrite the selection part of pmat with the new par_sel

p_index = (par_sel==1); % find the parameters that need to be fitted
loglikmax = -1 * transfer(parshat(p_index),pmat); % use transfer to obtain max likelihood!
Xhat = parshat(parnum); % best value of the profiled parameter

parmin  = pmat2(parnum,3); % use min and max of parameter as provided in script
parmax  = pmat2(parnum,4); % taken from pmat2, so these may be on log-scale

flag_better  = 0; % keep track of whether a better value is found. so we can stop the analysis if needed
for j = 1:2 % calculate from profile down and up from the ML estimate
    if j == 1
        % First calculating from max. lik. estimate downwards to zero ...
        logliklowest = 0;
        proffact     = -1; % going down ...
        
        if verbose == 1
            [figh,ft] = make_fig(1,1); % initialise a figure of correct size to keep track of the progress
            h1 = gca; % remember the current axis number for plotting!
            set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            hold on
            title(['Called from: ',scriptnm, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
            plot(h1,Xhat,0,'ko','MarkerFaceColor','w') % plot the max lik value as a circle
            drawnow % update the plot
            
            % if the parameter was on log-scale, mention it in the axis label
            if pmat(parnum,5) == 0
                xlabel(h1,['parameter ',names{parnum},' (log-scale)'],ft.name,ft.label)
            else
                xlabel(h1,['parameter ',names{parnum}],ft.name,ft.label)
            end
            ylabel(h1,'minus 2x log-likelihood ratio',ft.name,ft.label)
            
            drawnow
        end
    else
        % Now starting calculations from max. lik. estimate upwards ...
        proffact = 1; % going up ...
        % remember values from the first run (downwards) and clear 'em
        XcollD = Xcoll;
        loglikratD = loglikrat;
        clear Xcoll loglikrat;
    end
    
    p_try     = pmat; % temporary parameter vector
    Xtry      = Xhat; % initialise the parameter to try on the ML estimate
    logliktmp = 0;
    i         = 1;
    loglik_disp = -1; % make sure there is an output
    
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
                if Fstep * abs(Xtry_old) > 1.5 * Fstep_abs && Fstep > 1.5 * Fstep_min,
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
            % (bounds specified in the script file)
            Xtry = max(Xtry,parmin);
            Xtry = min(Xtry,parmax);
            
            p_try(parnum,1) = Xtry ; % put new profile value into the total try vector
            p_fit = p_try(p_index,1); % these are the ones that must be estimated now
            
            % Robustness of the optimisation routines is a real issue with
            % profiling, unfortunately. Here, I implement two methods: use
            % simulated annealling (probably less sensitive to local
            % minima), or a method with simplex sub-optimisations. The
            % latter seems to be best (but that will depend on the settings
            % of annealing).
            
            if subann == 1
                
                opt_ann.InitTemp  = 5; % much higher initial temperature
                opt_ann.Verbosity = 0; % controls output on screen (no output here)
                opt_ann.StopTemp  = 1e-2; % crude temperature for stopping as Simplex will finish it off
                [parshattemp,~] = anneal(@transfer,p_fit,opt_ann,p_try); % first rough estimation
            
            else
                
                % Now do fminsearch to get new best fit, with profile parameter fixed
                % Do the optimisation at least twice: once rough, and once a little more
                % precise. Hopefully, this will avoid local minima, but still
                % be rapid enough!
           
                % minimisation is done with the reduced parameter set, other parameters in par are kept
                % to the fixed value. This means that parshat will also contain the limited set only.
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
                
                if verbose == 1
                    plot(h1,Xcoll(i),loglikrat(i),'k.') % plot em on the fly!
                    drawnow
                end
                
                if brkprof == 1 && flag_better == 1 && logliktmp > logliklowest
                    % We have found a better value than the best one
                    % already (flag_better=1), but now it's getting worse
                    % again, so can break off the run.
                    par_better = packunpack(2,0,pmat_better); % return the better estimate in par_better
                    if verbose == 1
                        fprintf('  \n');
                        fprintf('Stopping the analysis: the profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,-1*loglikmax);
                        fprintf('Parameter values for new optimum \n');
                        fprintf('=================================================================================\n');
                        print_par(pmat_better); % write the better estimate to screen/diary in a formatted way
                        fprintf('=================================================================================\n');
                        fprintf('  \n');
                    end
                    diary off  % close results.out
                    return     % return to function calling us
                end
                
                if logliktmp < logliklowest % if we found a lower likelihood, remember it!
                    if logliktmp < -0.01 % ignore tiny improvements before breaking the analysis
                        flag_better  = 1; % mark that we have located a better optimum
                    end
                    
                    logliklowest = logliktmp;
                    p_disp       = p_try;
                    loglik_disp  = logliknew;
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
end
fclose(fid); % close the file profiles_newopt.out for writing

%% Calculating the 95% interval and print it on-screen

if logliklowest < 0 && flag_better == 1 % if the profiler found a better optimum, and we need to break
    par_better = packunpack(2,0,pmat_better); % return the better estimate in par_better
end

if verbose == 1
    diary ('results.out') % collect screen output in the diary "results.out"
    
    if logliklowest < 0 % if the profiler found a better optimum, display it!
        par_better = packunpack(2,0,pmat_better); % return the better estimate in par_better
        % this is done here again, because of the flag_better (a slightly
        % better output should not make auto_profiles.m do another
        % optimisation, but it does not hurt to put it on screen).
        
        loglik_disp = -1*loglikmax+logliklowest/2; % cannot use old loglik_disp as it might have been overwritten with a -1!
        
        fprintf('\n  The profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,-1*loglikmax);
        fprintf('Parameter values for new optimum \n');
        fprintf('=================================================================================\n');
        print_par(pmat_better); % write the better estimate to screen in a formatted way
        fprintf('=================================================================================\n');
        fprintf('  \n');
    end
end

XcollU     = Xcoll;
loglikratU = loglikrat;

if verbose == 1
    disp('  ')
    disp('95% confidence interval from the profile')
    disp('=================================================================================')
end

% this should not be needed, but sometimes we get an Inf in there ...
loglikratU = min(loglikratU,1000); 
loglikratD = min(loglikratD,1000); 

prof_coll = [XcollD' loglikratD';XcollU' loglikratU']; % construct entire profile
prof_coll = sortrows(prof_coll,1); % sort on basis of x-axis
% prof_coll will also be returned as optional output (log pars are on log-scale)

Xing = calc_xing(prof_coll,chicrit);  % find crossings of the chicrit
% back-transform Xing if needed for CIs on screen
if pmat(parnum,5) == 0
    Xing(:,1) = 10.^Xing(:,1);
end
if verbose == 1
    % column 2 is for warnings when hitting a boundary
    if Xing(1,2) == 1 % lowest crossing an upper boundary ...
        disp('Warning: taking the lowest parameter value as lower confidence limit')
    end
    if Xing(end,2) == 1 % highest crossing a lower boundary ...
        disp('Warning: taking the highest parameter value as upper confidence limit')
    end
    if size(Xing,1) == 2 % no problems, it is a single interval
        fprintf('%-10s interval: %#10.4g - %#1.4g \n',names{parnum},Xing(1,1),Xing(2,1))
    else
        disp('The confidence interval is a broken set (check likelihood profile to check these figures)')
        fprintf('%-10s interval: \n',names{parnum})
        for ix = 1:size(Xing,1)/2
            fprintf('   interval %1.1g: %#10.4g - %#1.4g \n',ix,Xing((ix-1)*2+1,1),Xing((ix-1)*2+2,1))
        end
    end
    
    % And make a nice plot of the profile
    plot(h1,prof_coll(:,1),prof_coll(:,2),'k-')
    
    % Plot the chi-square criterion for 95% and 1 degree of freedom
    minX = prof_coll(1,1); 
    maxX = prof_coll(end,1);
    plot(h1,[minX maxX],[chicrit chicrit],'k:')
    xlim(h1,[minX-0.05*(maxX-minX) maxX+0.05*(maxX-minX)]) % limit x-axis
    drawnow
    
    % set the y-axis to allow to read the graph when lik-ratio becomes very large
    plotmin = min(prof_coll(:,2));
    plotmax = max(10,min(30,max(prof_coll(:,2)))); % force plotmax between 10 and 30
    ylim(h1,[plotmin plotmax])
    
    plot(h1,Xhat,0,'ko','MarkerFaceColor','w') % plot the max lik value as a circle
        
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['profile_',scriptnm,'_par_',names{parnum}];%
        save_plot(gcf,savenm);
    end
    
    disp(' ')
    disp(['Time required: ' secs2hms(toc)])
    diary off  % close results.out
end

