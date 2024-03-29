function [ECx,ECx_lo,ECx_hi,Feff,ind_traits] = calc_ecx(par_plot,Tend,opt_ecx,opt_conf)

% Usage: [ECx,ECx_lo,ECx_hi,Feff,ind_traits] = calc_ecx(par_plot,Tend,opt_ecx,opt_conf)
% 
% Calculate ECx,t for all available traits with confidence intervals. This
% function should work with every TKTD model you throw at it, as long as
% there is at elast one state variables indicated with one of the dedicated
% traits: <glo.locS>, <glo.locL>, <glo.locR> (more may be added in the
% future).
% 
% Some calculation speed can be gained by adding an option that skips
% recalculation of the control response when running through a sample for
% CIs. At least, when only the tox parameters are fitted!
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%            structure from the saved sample is used
% <Tend>       time points (vector) at which to calculate ECx
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals
% 
% Author     : Tjalling Jager 
% Date       : February 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 X0mat

names     = glo2.names;
filenm    = glo.basenm;

X0mat_rem = X0mat; % remember X0mat before we change it
glo_rem   = glo;   % remember glo before we change it

backhaz   = opt_ecx.backhaz; % parameter name (as string) to set to zero for LCx/LPx calculation to remove background mortality
setzero   = opt_ecx.setzero; % parameter names (as string array) for extra paramaters to be set to zero
Feff      = opt_ecx.Feff;    % effect level (>0 en <1), x/100 in LCx (also used here for ECx)
ECx_plot  = opt_ecx.plot;    % set to 0 to NOT make a plot of LCx vs time
X_excl    = opt_ecx.statsup; % states to suppress from the calculations (e.g., locS)
par_read  = opt_ecx.par_read; % when set to 1 read parameters from saved set, but do NOT make CIs
notitle   = opt_ecx.notitle;  % set to 1 to suppress titles above plots

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
end

Tend = unique(Tend); % only unique ones and sort
t    = unique([0:max(Tend) Tend]); % time vector for running the model
if length(t) < 3
    t = unique([0 t(end)/2 t]); % just make sure there are at least 3
end
t = t(:); % make it a column
[~,loc_T] = ismember(Tend,t); % find where Tend values are in the total t vector
    
% see which other states are there
locS = [];
locL = [];
locR = [];
if isfield(glo,'locS') && ~ismember(glo.locS,X_excl) % then we have a state of survival
    locS = glo.locS; % collect the location
end
if isfield(glo,'locL') && ~ismember(glo.locL,X_excl) % then we have a state of body length
    locL = glo.locL; % collect the location
end
if isfield(glo,'locR') && ~ismember(glo.locR,X_excl) % then we have a state of reproduction
    locR = glo.locR; % collect the location
end

% And also add the states for the GUTS immobility package. For now, healthy
% only, since for that trait, it is easy to calculate ECx relative to the
% control (for death and immobile, the control is zero). There is a way to
% calculate ECx for death, but that would require summing healthy and
% immobile animals before calculating the effect (or take 1-death), so that
% is a bit more work.
loc_h = [];
if isfield(glo,'loc_h') % then we have a state of healthy
    loc_h = glo.loc_h; % collect the locations
end

ind_traits = [locS locL locR loc_h]; % indices for the traits we want from Xout
% We seem to have no interest for damage ...

% vector with initial values for the states in the simulations
X0mat_tmp    = X0mat(:,1); % take first column for our analysis
X0mat_tmp(1) = 0;          % start with 0 concentration

make_scen(-5,-1); % remove all spline info for exposure profiles as we only need constant exposure for ECx

% If we need CIs, load the best parameter set and the random sample from file
if type_conf > 0 || isempty(par_plot) % also if par_plot is not provided
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) == 1 % that means that no sample was found
        type_conf = -1; % no need to produce an error, just do analysis without CI
    end
    if isempty(par_plot) % if no par structure was entered in this function ...
        par_plot = par; % simply use the one from the sample file
    end
    if type_conf < 1 || par_read == 1 % then we don't want to make CIs
        type_conf = 0;  % don't make CIs anymore (when triggered by par_read)
        rnd       = []; % make sample empty
    end
end

% if ~isfield(par_plot,'tag_fitted') % apparently, parameters have not been fitted
%     warning('off','backtrace')
%     warning('You did not fit any parameters, so LCx or LPx is based on the values in the initial parameter matrix par.')
%     warning('Any CIs are made from the saved set in the MAT file.')
%     disp(' '), warning('on','backtrace')
% end

% backhaz is by default set to 'hb' in prelim_checks.m for the GUTS package
% identify background hazard and set it to zero
if isempty(backhaz) || ~isfield(par_plot,backhaz) % we need to make background mortality zero
    error('The function calc_ecx expects a parameter name in opt_ecx.backhaz, which matches a parameter name in your parameter structure that can be set zero to remove background mortality.')
else
    eval(['par_plot.',backhaz,'(1) = 0;']); % set parameter to zero in par_plot
    loc_zero = strcmp(names,backhaz)==1; % where is this parameter in the par structure?
end
% allow extra parameters to be set to zero, such as initial concentrations
if ~isempty(setzero)
    if ~iscell(setzero) % just to make sure it is a cell array
        setzero = {setzero}; % turn it into a cell array with one element
    end
    for i = 1:length(setzero)
        eval(['par_plot.',setzero{i},'(1) = 0;']); % set parameter to zero in par_plot
        loc_zero = loc_zero == 1 | strcmp(names,setzero{i})==1; % add this parameter to loc_zero
    end
end

%% Rough exploration of the concentration range
% To find out if there is an ECx,t, and where it approximately is. This
% should give us good starting ranges for all traits and all time points.

Xout  = call_deri(t,par_plot,X0mat_tmp); % use call_deri.m to provide the output for one scenario
Xctrl = Xout(loc_T,ind_traits);          % remember the relevant control output for the traits

c = 1; % start with concentration 1
Xout      = call_deri(t,par_plot,[c;X0mat_tmp(2:end)]); % use call_deri.m to provide the output for one scenario
Xout      = Xout(loc_T,ind_traits) ./ Xctrl; % only use the relevant relative output for the traits
Xout_coll = [c*ones(length(Tend),1) Tend(:) Xout]; % remember the relative output for the traits
% Xout_coll has two columns at the start for concentration and time
Xout1 = Xout; % remember the one at c=1

% while ~all(min(Xout_coll(:,3:end),[],2)<1-max(Feff)) && c < 1e6 % stop increasing c until there is large enough effect for all traits
% while ~all(min(Xout,[],1)<1-max(Feff)) && c < 1e6 % stop increasing c until there is large enough effect for all traits
while ~all(Xout(~isnan(Xout)) < 1-max(Feff)) && c < 1e6 % stop increasing c until there is large enough effect for all traits and all time points
    c = c * 10;
    Xout = call_deri(t,par_plot,[c;X0mat_tmp(2:end)]); % use call_deri.m to provide the output for one scenario
    Xout = Xout(loc_T,ind_traits) ./ Xctrl; % only use the relevant relative output for the traits
    Xout_coll = cat(1,Xout_coll,[c*ones(length(Tend),1) Tend(:) Xout]);
end

Xout = Xout1;
c = 1; % start again from concentration 1
% while ~all(max(Xout_coll(:,3:end),[],2)>1-min(Feff)) && c > 1e-7 % stop decreasing c until there is small enough effect for all traits
% while ~all(max(Xout(:,3:end),[],1)>1-min(Feff)) && c > 1e-7 % stop decreasing c until there is small enough effect for all traits
while ~all(Xout(~isnan(Xout)) > 1-min(Feff)) && c > 1e-7 % stop decreasing c until there is small enough effect for all traits and all time points
    c = c / 10;
    Xout = call_deri(t,par_plot,[c;X0mat_tmp(2:end)]); % use call_deri.m to provide the output for one scenario
    Xout = Xout(loc_T,ind_traits) ./ Xctrl; % only use the relevant relative output for the traits
    Xout_coll = cat(1,[c*ones(length(Tend),1) Tend(:) Xout],Xout_coll);
end

% see if this ranges catches all ECx,t
disp(' ')
remX = [];
for i_X = 1:length(ind_traits)
    if min(Xout_coll(:,2+i_X)) > 1-max(Feff) || max(Xout_coll(:,2+i_X)) < 1-min(Feff)
        if all(Xout_coll(:,2+i_X) > 1-min(Feff))
            % now only remove a state when, at no point at all, is there
            % enough effect for the smallest effect level. In other cases,
            % there may be at least enough to calculate some effect levels.
            remX = cat(2,remX,i_X); % remember that trait for removal
        end
        switch ind_traits(i_X)
            case locS
                disp('For survival, there is insufficient range of effects to calculate all ECx,t.')
            case locL
                disp('For body length, there is insufficient range of effects to calculate all ECx,t.')
            case locR
                disp('For reproduction, there is insufficient range of effects to calculate all ECx,t.')
            case loc_h
                disp('For healthy animals, there is insufficient range of effects to calculate all ECx,t.')
        end
    end
end

ind_traits(remX)    = []; % remove that trait from the trait list
Xout_coll(:,2+remX) = []; % remove that trait from the collected values
Xctrl(:,remX)       = []; % remove that trait from the control values

%% Calculate ECx,t exactly with fzero

f = waitbar(0,'Calculating ECx. Please wait.','Name','calc_ecx.m');

for i_T = 1:length(Tend) % run through time points
    Xout_tmp  = Xout_coll(Xout_coll(:,2)==Tend(i_T),:); % only keep the results for this time point
    Xctrl_tmp = Xctrl(i_T,:); % only keep controls for this time point
    
    for i_X = 1:length(ind_traits) % run through traits
        
        waitbar(((i_T-1)*length(ind_traits)+i_X)/(length(Tend)*length(ind_traits)),f); % update waiting bar
        
        for i_F = 1:length(Feff) % run through effect levels
            
            ind_2    = find(Xout_tmp(:,i_X+2)<1-Feff(i_F),1,'first');
            ind_1    = find(Xout_tmp(:,i_X+2)>1-Feff(i_F),1,'last');
            EC_range = Xout_tmp([ind_1 ind_2],1); % range where ECx,t is located
            if numel(EC_range) == 2
                % use fzero to zero in on the exact value
                ECx{i_F}(i_T,i_X) = fzero(@calc_ecx_sub,EC_range,[],Tend(i_T),par_plot,Feff(i_F),X0mat_tmp(2:end),Xctrl_tmp(i_X),ind_traits(i_X)); % find the ECx,t
            else
                ECx{i_F}(i_T,i_X) = NaN; % then a proper range was not found
            end
        end
    end
end

close(f) % close the waiting bar

%% Display results without CI on screen

if type_conf > 0 % only do this when we'll go into CI calculation next
    
    % disp(' ')
    disp('Results for ECx without CIs (they are calculated next)')
    disp('================================================================================')
    fprintf('%s    ',glo.xlab)
    for i_F = 1:length(Feff)
        fprintf('EC%1.0f    ',100*Feff(i_F))
    end
    fprintf('(%s)',glo.leglab2)
    fprintf('\n')
    disp('================================================================================')
    
    for i_X = 1:length(ind_traits) % run through traits
        switch ind_traits(i_X)
            case locS
                disp('  Survival    :')
            case locL
                disp('  Body length :')
            case locR
                disp('  Reproduction:')
            case loc_h
                disp('  Healthy     :')
        end
        
        for i_T = 1:length(Tend) % run through time points
            fprintf('%3.0f ',Tend(i_T))
            for i_F = 1:length(Feff) % run through effect levels
                fprintf('%#10.3g     ',ECx{i_F}(i_T,i_X))
            end
            fprintf('\n')
        end
        disp('================================================================================')
    end
end

%% Calculate confidence intervals on the ECx,t

% initialise matrices to catch highest and lowest results from
% sample, per state and per time point
for i_F = 1:length(Feff)
    ECx_lo{i_F} = nan(length(Tend),length(ind_traits));
    ECx_hi{i_F} = nan(length(Tend),length(ind_traits));
end

if type_conf > 0 % if we make CIs ...
    n_sets   = size(rnd,1); % number of samples from parameter space
    pmat     = packunpack(1,par,0); % transform structure *from saved set* into a regular matrix
    % it is better to use the saved par, as there may be differences in the
    % log-setting of parameters between the saved set and the optimised
    % par_out matrix (especially when using the alllog option in
    % calc_slice).
    
    par_comp(par,par_plot,cat(2,backhaz,setzero)) % compare par from input with the one from the MAT file
    ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
    ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!
    
    f = waitbar(0,'Calculating confidence intervals on ECx. Please wait.','Name','calc_ecx.m');
    
    for k = 1:n_sets % run through all sets in the sample
        
        waitbar(k/n_sets,f); % update waiting bar
        
        pmat(ind_fit,1) = rnd(k,:); % replace values in pmat with the k-th random sample from the MCMC
        % put parameters that need to be fitted on log scale back on normal
        % scale (as call_deri requires normal scale, in contrast to transfer.m)
        if sum(ind_logfit)>0
            pmat(ind_logfit,1) = 10.^(pmat(ind_logfit,1));
        end
        % Note: pmat is on normal scale here, but the sample in rnd contains
        % the value on a log scale, if a parameter is fitted on log scale.
        pmat(loc_zero,1) = 0; % make parameter for background mortality (and possibly others) zero in each set of the sample!
        % do this after the back-transformation step.
        par_k = packunpack(2,0,pmat); % transform parameter matrix into a structure
    
        % Calculate the control response for this parameter set. When only
        % the tox parameters are fitted, this is superfluous. However, we
        % should not make a priori assumptions about how this function will
        % be used! E.g., for GUTS cases, hb may be fitted as well, and we
        % need to have the effect relative to the control for THIS set of
        % the sample.
        Xout  = call_deri(t,par_k,X0mat_tmp); % use call_deri.m to provide the output for one scenario
        Xctrl = Xout(loc_T,ind_traits);       % remember the relevant control output for the traits
        
        for i_T = 1:length(Tend) % run through time points
            
            for i_X = 1:length(ind_traits) % run through traits
                Xctrl_tmp = Xctrl(i_T,i_X); % only keep controls for this time point and trait
                
                for i_F = 1:length(Feff) % run through effect levels
                    if ~isnan(ECx{i_F}(i_T,i_X))
                        ECx_tmp = fzero(@calc_ecx_sub,ECx{i_F}(i_T,i_X),[],Tend(i_T),par_k,Feff(i_F),X0mat_tmp(2:end),Xctrl_tmp,ind_traits(i_X)); % find the ECx,t
                        % collect the min and max found so far
                        ECx_lo{i_F}(i_T,i_X) = min(ECx_lo{i_F}(i_T,i_X),ECx_tmp);
                        ECx_hi{i_F}(i_T,i_X) = max(ECx_hi{i_F}(i_T,i_X),ECx_tmp);
                    end
                end
            end
        end
    end
    
    close(f) % close the waiting bar

end

diary('results.out') % collect output in the diary "results.out"
% disp(' ')
switch type_conf
    case 0
        disp('ECx without confidence intervals');
    case 1
        disp('ECx with CIs: Bayesian 95% credible interval');
    case 2
        disp('ECx with CIs: 95% pred. likelihood, shooting method');
    case 3
        disp('ECx with CIs: 95% pred. likelihood, parspace explorer');
end
disp('================================================================================')
fprintf('%s    ',glo.xlab)
for i_F = 1:length(Feff)
    fprintf('EC%1.0f (95 perc. CI)    ',100*Feff(i_F))
end
fprintf('(%s)',glo.leglab2)
fprintf('\n')
disp('================================================================================')

for i_X = 1:length(ind_traits) % run through traits
    ECx_max(i_X) = NaN;
    
    switch ind_traits(i_X)
        case locS
            disp('  Survival    :')
        case locL
            disp('  Body length :')
        case locR
            disp('  Reproduction:')
        case loc_h
            disp('  Healthy     :')
    end
    
    for i_T = 1:length(Tend) % run through time points
        fprintf('%3.0f ',Tend(i_T))
        for i_F = 1:length(Feff) % run through effect levels
            fprintf('%#10.3g (%#10.3g - %#10.3g) ',ECx{i_F}(i_T,i_X),ECx_lo{i_F}(i_T,i_X),ECx_hi{i_F}(i_T,i_X))
            ECx_max(i_X) = max([ECx_max(i_X),ECx{i_F}(i_T,i_X),ECx_hi{i_F}(i_T,i_X)]);
        end
        fprintf('\n')
    end
    disp('================================================================================')

end
diary off

%% Make a plot of ECx versus time

if ECx_plot ~= 0
    
    n = length(Feff);
    m = length(ind_traits);
    [~,ft] = make_fig(m,n,2); % create a figure window of correct size
    
    for i_X = 1:length(ind_traits) % run through traits
        
        for i_F = 1:length(Feff) % run through effect levels
            
            h_pl = subplot(m,n,(i_X-1)*length(Feff)+i_F); % make a sub-plot
            hold on
            h_ax = gca;
            set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
                p = get(h_pl,'position');
                p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
                set(h_pl, 'position', p);
            end
            
            % create axis labels
            if i_X == length(ind_traits) % only x-label at last row
                xlab = [glo.xlab]; % label
                xlabel(xlab,ft.name,ft.label)
            else
                set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
            end
            if i_F == 1 % only y-label in first column
                ylab = ['ECx ',glo.ylab{ind_traits(i_X)}]; % label
                ylabel(ylab,ft.name,ft.label)
            else
                set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
            end
            if i_X == 1 % only title in first row
                title(['EC',num2str(100*Feff(i_F)),' (',glo.leglab2,')'])
            end
            
            if type_conf > 0 % then we have CIs to plot
                % Little trick to fill the area between the two curves, to
                % obtain a coloured confidence interval as a band.
                ind_ok = ~isnan(ECx{i_F}(:,i_X)); % there may be NaNs in there, that we don't use
                t2  = [Tend(ind_ok)';flipud(Tend(ind_ok)')]; % make a new time vector that is old one, plus flipped one
                Xin = [ECx_lo{i_F}(ind_ok,i_X);flipud(ECx_hi{i_F}(ind_ok,i_X))]; % do the same for the plot line, hi and lo
                fill(t2,Xin,'g','LineStyle','none','FaceAlpha',1) % and fill this object
            end
            
            plot(Tend,ECx{i_F}(:,i_X),'ko-','MarkerFaceColor','y') % plot ECx values
            ylim([0 1.05*ECx_max(i_X)]) % limit y-axis
            xlim([0 Tend(end)])
        end
    end
    
    if notitle == 0
        switch type_conf
            case 0
                tit = 'No confidence intervals';
            case 1
                tit = 'CIs: Bayesian 95% credible interval';
            case 2
                tit = 'CIs: 95% pred. likelihood, shooting method';
            case 3
                tit = 'CIs: 95% pred. likelihood, parspace explorer';
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,tit,'HorizontalAlignment','center','VerticalAlignment', 'top');
    end
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['epx_window_plot_',filenm];%
        save_plot(figh,savenm);
    end
    
end

% return the globals to their original states
X0mat = X0mat_rem;
glo   = glo_rem;

function crit = calc_ecx_sub(c,Tend,par_out,Feff,X0,Xctrl,locX)

% Usage: crit = calc_ecx_sub(c,Tend,par_out,Feff,X0,Xctrl,locX)
%
% This function calculates the criterion to be used by <fzero> to calculate
% the ECx. This function is generally useful for TKTD models, including
% GUTS and DEBtox. 
%
% Inputs:
% <c>       the concentration to try to see if it is the ECx
% <Tend>    the time at which to calculate the ECx
% <par_out> the optimised parameter set
% <Feff>    the fraction effect (x/100 in ECx)
% <X0>      initial values for the state variables
% <Xctrl>   value for the control trait (as ECx is relative to control)
% <locX>    location of trait of interest in Xout

Xtst2 = call_deri([0 Tend/2 Tend],par_out,[c;X0]); % response at conc. c
crit  = (Xtst2(end,locX)/Xctrl)-(1-Feff); % zero when end value is x*100% effect



