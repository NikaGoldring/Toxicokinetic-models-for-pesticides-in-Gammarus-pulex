function plot_tktd(par,opt_tktd,opt_conf)

% Usage: plot_tktd(par,opt_tktd,opt_conf)
%
% This function provides a series of additional plots, specifically for
% TKTD purposes (analysis of toxic effects data over time). The standard
% plotting routine always plots traits versus time, whereas other
% plotting options may be more useful: instead of plotting all data in a
% single plot, this script provides multiplots (different plots for
% different treatments).
% 
% At this moment, this function is used for GUTS and DEBtox analysis,
% including analysis of binary mixtures (although that does require some
% care and further checking), and including the future immobility package.
%
% Note that the error bars plotted on observations are approximate 95% CIs.
% For survival, this is the Wilson score and for continuous traits 2xSE.
% Leave opt_conf empty to plot without CIs on model curves. If par is left
% empty, the function attempts to load it from a saved sample.
%
% Author     : Tjalling Jager
% Date       : September 2020
% Web support: http://www.debtox.info/byom.html
%
% LIMITATIONS
% - Works with removed/missing animals (for survival), but this needs more
% checking to make sure it works in every case
% - Sampling error is not used yet.
% - Works for binary mixtures, but this also needs more checking

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global DATA W X0mat % make the data set and initial states global variables
global glo glo2     % allow for global parameters in structure glo

plot_repls   = opt_tktd.repls;   % plot individual replicates (1) or means (0)
plot_obspred = opt_tktd.obspred; % plot predicted-observed plots (1) or not (0)
plot_preds   = opt_tktd.preds;   % set to 1 only plot predictions from X0mat without data
plot_addzero = opt_tktd.addzero; % set to 1 to always add a concentration zero to X0mat
max_on_C     = opt_tktd.max_exp; % set to 1 to maximise exposure/damage plots on exposure rather than damage
notitle      = opt_tktd.notitle; % set to 1 to suppress titles above plots
transf       = opt_tktd.transf;  % set to 1 to calculate means and SEs including transformations
plot_min     = opt_tktd.min;     % set to 1 to show a dotted line for the control (lowest) treatment
statsup      = opt_tktd.statsup; % states to suppress from the plots (e.g., locS)

% some hidden options that may be useful in some cases
lim_data   = 0; % set to 1 to limit y-axis to data (default)
rem_y_info = 0; % set to 1 to remove all y-axis labels from first column (only when plotting data for now)

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
    opt_conf.sens = 0; % no sensitivities when calling calc_conf
end

n_D     = glo2.n_D; % number of data sets per state
mod_t   = glo.t;    % time vector for model curves
mod_t   = mod_t(:); % make sure it is a column vector
leglab1 = glo.leglab1; % legend label 1
leglab2 = glo.leglab2; % legend label 2

% We need to take care of slow kinetics, but I have used two formats for
% the global glo.fastslow: I used to use a character in there, but I like
% to move to numbers, and in fact a 2-element vector (to also allow for
% fast damage repair).
slowkin = 0; % flag for special case of slow kinetics!
if isfield(glo,'fastslow')
    if ischar(glo.fastslow) && glo.fastslow == 's'
        slowkin = 1;
    elseif glo.fastslow(1) == 2
        slowkin = 1;
    end
end

% incl_conf = 0; % do NOT include CIs on plots, unless ...
% % incl_samerr = 0; % do NOT include sampling error bounds
% if ~isempty(opt_conf) % if left empty, you want to skip all CIs
%     incl_conf = 1; % include CIs on plots (but there must be a sample saved)
%     opt_conf.sens = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
% %     if opt_conf.samerr == 1 % if you calculated sampling errors before, I assume you want them here too ...
% %         incl_samerr = 1; % include sampling error bounds
% %     end
% end 

if isempty(par) % if par is not provided, load it from the saved sample (if present)
    if ~isempty(opt_conf)
        [~,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    else
        error('No parameter structure was entered, and also no info in opt_conf, so cannot continue.')
    end
elseif ~isempty(opt_conf) && ~isempty(opt_conf.set_zero)
    % we may want to make a parameter zero (esp. background mortality)
    % eval(['par.',opt_conf.set_zero,'(1) = 0;']); % set parameter to zero in par

    setzero = opt_conf.set_zero; % this is a string or cell array of strings with parameter names
    if ~iscell(setzero) % just to make sure it will be a cell
        setzero = {setzero}; % turn it into a cell array with 1 element
    end
    for i = 1:length(setzero)
        eval(['par.',setzero{i},'(1) = 0;']); % set parameter to zero in par
    end
    
    % This is a bit annoying ... load_rnd will set hb to zero, but only if
    % the par is read from the saved file. If it is entered as input to
    % this function, we need to do it here ... For the sample, it is
    % automatically done by calc_conf.
end

mod_c = X0mat(1,:); % all scenarios that are to be plotted
if ~ismember(0,mod_c) && plot_addzero == 1 % make sure that there is a zero in there!
    % this does slow down calculations, especially with CIs ...
    mod_c = cat(2,0,mod_c); % add a zero
    X0mat = cat(2,[0;X0mat(2:end,1)],X0mat); % also add the control to X0mat
    % NOTE: this assumes that control will have same initial values for the
    % states as the first entry in X0mat!
end

% Deaths per interval (both observed and predicted) are collected, but not
% used yet. For the immobility case, this will be nonsense.
%
% Multiple data sets for the exposure is now implemented, for BINARY
% mixtures only! (and assuming coding identifiers with a factor of 10).
% This still needs serious testing.

%% Collect locations of trait states in the state-variable list

% See which other states are there
locD = [];
locS = [];
locL = [];
locR = [];
locC = [];
if isfield(glo,'locS') % then we have a state of survival
    locS = glo.locS; % collect the locations
end
if isfield(glo,'locL') % then we have a state of body length
    locL = glo.locL; % collect the locations
end
if isfield(glo,'locR') % then we have a state of reproduction
    locR = glo.locR; % collect the locations
end
if isfield(glo,'locC') % then we have a state of internal concentration
    locC = glo.locC; % collect the locations
end
if isfield(glo,'locD') % then we have a state of damage
    locD = glo.locD; % collect the locations
end
% There can be multiple damage states (e.g., for mixtures)
if length(locD) > 2
    error('At this moment, the plotting routine in plot_tktd.m can not accommodate more than two compounds/damage levels.')
end

% And also add the states for the GUTS immobility package.
loc_h = [];
if isfield(glo,'loc_h') % then we have a state of healthy
    loc_h = glo.loc_h; % collect the locations
end
loc_i = [];
if isfield(glo,'loc_i') % then we have a state of immobile
    loc_i = glo.loc_i; % collect the locations
end
loc_d = [];
if isfield(glo,'loc_d') % then we have a state of dead
    loc_d = glo.loc_d; % collect the locations
end
loc_id = [];
if isfield(glo,'loc_id') % then we have a state of immobile or dead
    loc_id = glo.loc_id; % collect the locations
end

%% Initialise matrices to collect model output

for i_D = 1:length(locD)
    XoutD{i_D} = nan(length(mod_t),length(mod_c)); % matrices to collect damage model curves
end
% Now assume that each damage level also has a separate exposure profile!
for i_D = 1:length(locD)
    XoutC{i_D} = nan(length(mod_t),length(mod_c)); % matrix to collect external concentrations
end
locP  = [locS loc_h loc_i loc_d loc_id]; % collect output states for quantal states
locPs = [loc_i loc_d loc_id]; % collect output states for quantal states that need some special attention
locX  = [locC locP locL locR]; % collects output states for all traits
% NOTE: I just add the internal concentration to the states. That is not
% really nice, but we won't usually have that state in TKTD models.

% at this point, remove the states that we were asked to suppress
[a,b]  = ismember(statsup,locP);
b(b==0) = [];
locP(b) = [];
[a,b]  = ismember(statsup,locPs);
b(b==0) = [];
locPs(b) = [];
[a,b]  = ismember(statsup,locX);
b(b==0) = [];
locX(b) = [];
[a,b]  = ismember(statsup,locD);
b(b==0) = [];
locD(b) = [];

% There can be multiple output states
for i_X = 1:length(locX)
    XoutX{i_X} = nan(length(mod_t),length(mod_c)); % matrices to collect model curves for other states
end

%% Calculate and collect model output

modmaxX = zeros(1,length(locX)); % catch overall maximum model curves on y-axis
modmaxD = zeros(1,length(locD)); % catch overall maximum model curves on y-axis

for i_c = 1:length(mod_c) % run through modelled scenarios in X0mat
    
    Xout = call_deri(mod_t,par,X0mat(:,i_c)); % use call_deri.m to provide the model output for each scenario
     
    for i_X = 1:length(locX) % run through other states
        XoutX{i_X}(:,i_c) = Xout(:,locX(i_X)); % store model curve in correct place
        modmaxX(i_X) = max(modmaxX(i_X),max(XoutX{i_X}(:,i_c))); % update maximum of states
    end
    
    for i_D = 1:length(locD) % run through damage states
        XoutD{i_D}(:,i_c) = Xout(:,locD(i_D)); % store model curve in correct place
        modmaxD(i_D) = max(modmaxD(i_D),max(XoutD{i_D}(:,i_c))); % update maximum of damage
    end
    
    for i_D = 1:length(locD) % run through damage states (assume that each damage state has an exposure profile)
        XoutC{i_D}(:,i_c) = mod_c(i_c) * ones(length(mod_t),1); % if no exposure profile is specified, simply copy c over all time points
        
        if isfield(glo,'int_scen') % if int_scen exists: use it to derive current external conc.
            
            if length(locD) == 1 % in regular cases, there is only 1 damage
                if ismember(mod_c(i_c),glo.int_scen) % is this c in the scenario range global?
                    XoutC{i_D}(:,i_c) = make_scen(-3,mod_c(i_c),mod_t); % -3 to ask for concentration vector
                end
            else % assume it is a BINARY mixture
                % Extract identifiers for each compound (assume a factor of 10 was used in coding the identifiers)
                cI(1) = 10*floor(mod_c(i_c)/10); % identifier for compound A
                cI(2) = mod_c(i_c) - cI(1);      % identifier for compound B
                
                if ismember(cI(i_D),glo.int_scen) % is cI in the scenario range global?
                    XoutC{i_D}(:,i_c) = make_scen(-3,cI(i_D),mod_t); % -3 to ask for concentration vector
                else
                    XoutC{i_D}(:,i_c) = cI(i_D) * ones(length(mod_t),1); % if no exposure profile is specified, simply copy c over all time points
                    % this will mainly be needed for the single exposures,
                    % where the other chemical is zero.
                end
            end
            
        end
    end
end

% Initialise matrices to catch confidence intervals
for i_X = 1:length(locX) % run through other states
    XloX{i_X} = nan(length(mod_t),length(mod_c));
    XhiX{i_X} = nan(length(mod_t),length(mod_c));
end
for i_D = 1:length(locD) % run through damage states
    XloD{i_D} = nan(length(mod_t),length(mod_c));
    XhiD{i_D} = nan(length(mod_t),length(mod_c));
end

% Calculate CIs on the model curves and collect them in the correct
% matrices. The function calc_conf is used to calculate CIs.
if type_conf > 0
    out_conf = calc_conf(par,opt_conf,1); % TEMP ALWAYS PLOT REPLICATES (should not matter as we don't plot sampling error yet)
    for i_X = 1:length(locX)
        XloX{i_X} = out_conf{1}{locX(i_X)};
        XhiX{i_X} = out_conf{2}{locX(i_X)};
        modmaxX(i_X) = max(modmaxX(i_X),max(XhiX{i_X}(:))); % update maximum of states
    end
    for i_D = 1:length(locD)
        XloD{i_D} = out_conf{1}{locD(i_D)};
        XhiD{i_D} = out_conf{2}{locD(i_D)};
        modmaxD(i_D) = max(modmaxD(i_D),max(XhiD{i_D}(:))); % update maximum of states
    end
end

for i_X = 1:length(locX)
    if ismember(locX(i_X),locP) % if it is a probability state ...
        modmaxX(i_X) = 1; % maximum can always be one
    end
end
    
%% Plot predictions without data
% This simply plots the model curves for all treatments specified in X0mat.
% This is handy as this function focusses strongly on data: the data are
% leading, and only scenarios are plotted that have data. This section uses
% the option plot_preds to bypass that and plot all scenarios in X0mat, but
% without any data. This section of the code is copied from the main
% plotting code below.

if plot_preds == 1
    
    % assume that the lowest c is the control ... as control will be plotted in all panels
    [~,loc_min] = min(mod_c); % where is the scenario with the minimum overall in model set?
    loc_min = loc_min(1); % just make sure there's only 1 (superfluous)
    
    m      = length(locD) + length(locX); % number of states to plot (number of rows in multiplot)
    n      = length(mod_c); % number of treatments to plot (columns in multiplot)
    [figh,ft] = make_fig(m,n,2); % create a figure window of correct size
    xlab   = glo.xlab; % text on x-axis
    
    for i_c = 1:n % run through relevant treatments in X0mat
        row = -1;
        
        % First plot damage and exposure scenario
        for i_D = 1:length(locD) % run through damage states (there may be more than one)
            row = row + 1; % go to next row of plots for damage
            
            if slowkin == 1 % for slow kinetics
                ylab = 'conc./damage time';
            else            
                ylab = glo.ylab{locD(i_D)}; % text on y axis
            end
            
            mod_plt = [mod_t XoutD{i_D}(:,i_c) XloD{i_D}(:,i_c) XhiD{i_D}(:,i_c)]; % model things to plot: mean, low CI and high CI
            h_pl(row*n+i_c) = plot_helper(m,n,row*n+i_c,ft,xlab,ylab,mod_plt,{[]}); % use helper function to plot model curves
            
            if opt_tktd.plotexp == 1 % if we want to plot the exposure ...
                area(mod_t,XoutC{i_D}(:,i_c),'FaceColor','b','FaceAlpha',0.1); % plot concentration profile as a filled area
            end
            
            if i_D == 1 && opt_tktd.flip == 0 % only place title in first row
                % create title above panel (but not if we intend to flip later!)
                if n_D>1 % if there is more than 1 data set
                    tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
                else
                    tit = '';
                end
                if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
                    title([tit,Ltmp],ft.name,ft.title)
                else
                    title([tit,leglab1,num2str(mod_c(i_c)),' ',leglab2],ft.name,ft.title)
                end
            end
            xlim([0 mod_t(end)]) % limit x-axis to model time vector
            
            if max_on_C == 1
                ylim([0 1.05*max(XoutC{i_D}(:))]) % limit y-axis to exposure concentration
            else
                ylim([0 1.05*modmaxD(i_D)]) % limit y-axis to calculated damage
            end
            
        end
        
        % Plot output states for traits
        for i_X = 1:length(locX) % run through all other state variables
            
            row = row + 1; % go to next row to plot next state variable
            % collect info for plotting model curve with CI (and plotting the control as well)
            if plot_min == 1
                mod_plt = [mod_t XoutX{i_X}(:,i_c) XloX{i_X}(:,i_c) XhiX{i_X}(:,i_c) XoutX{i_X}(:,loc_min)];
            else
                mod_plt = [mod_t XoutX{i_X}(:,i_c) XloX{i_X}(:,i_c) XhiX{i_X}(:,i_c)];
            end
            
            ylab    = glo.ylab{locX(i_X)}; % label for y-axis
            
            h_pl(row*n+i_c) = plot_helper(m,n,row*n+i_c,ft,xlab,ylab,mod_plt,{[]}); % use helper function to plot model curves
            xlim([0 mod_t(end)]) % limit x-axis to model time vector
            ylim([0 1.05*max(modmaxX(i_X),1e-6)]) % limit y-axis to calculated state (avoid error when a max is zero)
            
            % TEST also plot exposure scenario for internal concentrations
            if ismember(locX(i_X),locC) && opt_tktd.plotexp == 1 % then we have internal concentrations to plot, and asked to plot exposure ...
                area(mod_t,XoutC{i_D}(:,i_c),'FaceColor','b','FaceAlpha',0.1); % plot concentration profile as a filled area (option 'FaceAlpha',0.25 removed)
                ylim([0 max(1e-6,max(1.05*modmaxX(i_X),1.05*max(XoutC{i_D}(:))))]) % limit y-axis to exposure level or concentration
            end
            
            % There may be cases where we want to use this function,
            % but without damage state. To make sure that titles are
            % plotted in the first row, we need to repeat this code,
            % which is a bit of a pain.
            if isempty(locD) || opt_tktd.flip == 1
                if i_X == 1 % only place title in first row
                    % create title above panel
                    if n_D>1 % if there is more than 1 data set
                        tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
                    else
                        tit = '';
                    end
                    if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                        Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
                        title([tit,Ltmp],ft.name,ft.title)
                    else
                        title([tit,leglab1,num2str(mod_c(i_c)),' ',leglab2],ft.name,ft.title)
                    end
                end
            end
            
        end
    end
    
    % Experimental! Try to swap the first and second row of the plot. This is
    % useful if the plot includes damage AND internal concentrations. Since
    % internal concentrations are part of locX they are plotted AFTER damage,
    % while BEFORE would seem more logical. Code below works, but the titles
    % above the columns now also shift one row down. That is why this option is
    % also used to decide when to plot the titles above. This is quite
    % cumbersome, but other solutions would require more structural revision of
    % the code.
    if opt_tktd.flip == 1
        hFig    = gcf;
        hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
        hAxes   = flipud(hAxes); % reverse order
        for i = 1:n % run through rows of plot
            h_tst = h_pl([n+i i]); % handle for plot on first row and second row, reversed!
            set(hAxes([(i-1)*m+1 (i-1)*m+2]),{'position'},{h_tst.Position}.') % use them to set position in hAxes
            % Note: hAxes has a different order for the plots than h_pl: it
            % goes through the rows rather than through the columns.
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
        savenm = ['tktd_preds_plot_',glo.basenm];%
        save_plot(figh,savenm);
    end
    
    return % rest of the code is for plotting data, so we can return
end

%% Break up data set keeping all treatments in the data set(s)
% This allows for plotting replicates as individual data points.

data.maxt = zeros(n_D,1); % catch overall maximum on x-axis
data.max  = zeros(length(locX),n_D); % catch overall maximum on y-axis

for i_X = 1:length(locX) % run through trait states for other states (not damage since there is no data anyway)
    
    for i_d = 1:n_D % run through all data sets per state
        
        plotdata = DATA{i_d,locX(i_X)}; % extract data matrix
        w        = W{i_d,locX(i_X)};    % extract weight factors
        
        data_tmp = plotdata(2:end,2:end); % extract data set for this state (without t column and c row)
        data.c{i_X,i_d} = plotdata(1,2:end); % conc vector in the data set i_d
        t_id = plotdata(2:end,1); % time vector in the data set i_d
        lam  = plotdata(1,1); % indicator for transformations
        
        % remove from data.c and data_tmp the treatments that are not in mod_c
        ind_ok = ismember(data.c{i_X,i_d},mod_c); 
        data.c{i_X,i_d} = data.c{i_X,i_d}(:,ind_ok);
        data_tmp = data_tmp(:,ind_ok);
        
        for i_c = 1:length(data.c{i_X,i_d}) % run through treatments in this data set
            
            data_tmp_ic = data_tmp(:,i_c); % data for this treatment only
            w_ic = w(:,i_c); % vector with missing animals for this treatment only
            
            if ismember(locX(i_X),locP) % then we have quantal data
                
                % Combine all vectors into the overall object, as a
                % nice matrix for data set <i_d> and treatment <i_c>. columns:
                % time, survivors, nr. deaths, surv. prob., CI min, CI max
                [data_out,S0] = recalc_data(data_tmp_ic,w_ic,lam);
                
                if ismember(locX(i_X),locPs) % than we have immobility and death states separately
                    data_out(1,:) = 0; % make all initial things at t=0 zero!
                    % Needed as first row in data for immobility and death
                    % states will be the total amount of animals with which
                    % the experiment started!
                end
                
                data.X{i_X,i_d,i_c} = [t_id data_out];
                
                % also collect the interpolated model values belonging to the data points!
                loc_c = find(mod_c == data.c{i_X,i_d}(i_c)); % where is this concentration in the total model vector from X0mat?
                data.mod{i_X,i_d,i_c} = [t_id interp1(mod_t,XoutX{i_X}(:,loc_c),t_id) interp1(mod_t,XloX{i_X}(:,loc_c),t_id) interp1(mod_t,XhiX{i_X}(:,loc_c),t_id)];
                
                % also collect the predicted deaths from the model survival
                % and the initial number of animals S0
                Dmod = data.mod{i_X,i_d,i_c}(:,[2:4]); % extract survival probabilties with CI
                Dmod = Dmod - [Dmod(2:end,:);zeros(1,3)]; % death probabilities in each interval, from the data
                data.moddeaths{i_X,i_d,i_c} = [t_id S0.*Dmod(:,1) S0 .*Dmod(:,2) S0.*Dmod(:,3)];
                
                data.max(i_X,i_d) = 1; % max(data.max(i_X,i_d),max(min(1,a+b))); % update overal maximum using highest edge of CI
                % probably also need to account for model value, which may be
                % higher (for other states than survival).
                data.maxt(i_d) = max(data.maxt(i_d),max(t_id)); % update overal maximum using highest edge of CI  
                
            else % then we have length or repro data ...
                
                w_ic(isnan(w_ic)) = 0; % replace any NaN weights by zero (this is for plotting outliers)
                
                % data_tmp_i_c is the vector with traits for this treatment only
                data.X{i_X,i_d,i_c} = [t_id data_tmp_ic nan(length(t_id),2) w_ic]; % store the time points and the state values
                % the NaNs are for the CI, which we don't have
                
                % also collect the interpolated model values belonging to the data points!
                loc_c = find(mod_c == data.c{i_X,i_d}(i_c)); % where is this concentration in the total model vector from X0mat?
                data.mod{i_X,i_d,i_c} = [t_id interp1(mod_t,XoutX{i_X}(:,loc_c),t_id) interp1(mod_t,XloX{i_X}(:,loc_c),t_id) interp1(mod_t,XhiX{i_X}(:,loc_c),t_id)];
                
                data.max(i_X,i_d) = max(data.max(i_X,i_d),max(data_tmp_ic)); % update overal maximum using highest value of state    
                data.maxt(i_d)    = max(data.maxt(i_d),max(t_id)); % update overal maximum using highest value of time vector
            end
        end
    end
end

%% Reconstruct data set as means!
% For survival, we can simply sum all replicates (taking care of possible
% missing/removed animals. When there are NaNs in one or more of the
% replicates, the mean will be NaN as well! There may be a workaround, but
% I don't see a simple solution at this moment.
% 
% For body length and reproduction, the mean is taken and standard error
% calculated. When a transformation is used in fitting, this may be a bit
% misleading.

means.maxt = zeros(n_D,1); % catch overall maximum on x-axis
means.max  = zeros(length(locX),n_D); % catch overall maximum on y-axis
flag_meanskip = 0; % flag if there is a mean skipped because there are NaNs in some replicates

for i_X = 1:length(locX) % run through trait states (not damage)
    
    for i_d = 1:n_D % run through all data sets per state
        
        [isinset] = ismember(mod_c,data.c{i_X,i_d}); % are the scenarios for this data set in mod_c?
        mod_c_i_d = mod_c(isinset); % only that part of mod_c that is used for data set i_d
        
        plotdata = DATA{i_d,locX(i_X)}; % extract data
        w        = W{i_d,locX(i_X)};    % extract weight factors
        
        data_tmp = plotdata(2:end,2:end); % extract data set for survivors
        t_id  = plotdata(2:end,1); % time vector in the data set i_d
        means.c{i_X,i_d} = mod_c_i_d; % unique concentration vector in the data set i_d
        c_tmp = plotdata(1,2:end); % conc vector in the data set i_d
        lam   = plotdata(1,1); % indicator for transformations
        
%         means.max(i_X,i_d) = 1; % make sure it is defined even when there are no data to plot
        

        for i_c = 1:length(mod_c_i_d) % run through unique treatments in data set
            
            [~,ind_i_c]  = ismember(c_tmp,mod_c_i_d(i_c)); % where in data set is this unique c? (can be more than 1!)
            data_tmp_ic = data_tmp(:,ind_i_c==1); % extract all replicates for this treatment
            w_ic = w(:,ind_i_c==1); % extract all weights for this treatment
            
            if ismember(locX(i_X),locP) % then we have quantal data
                
                % Combine all vectors into the overall object, as a
                % nice matrix for data set <i_d> and treatment <i_c>. columns:
                % time, survivors, nr. deaths, surv. prob., CI min, CI max
                [data_out,S0] = recalc_data(data_tmp_ic,w_ic,lam);
                
                chck_nans = mean(isnan(data_tmp_ic),2); % this is 0 when there are no NaNs, and 1 when there are NaNs in all replicates
                if ~all(chck_nans == 1 | chck_nans == 0) % if they are NOT all 0 or all 1, set a flag as we cannot plot some means!
                    flag_meanskip = 1; % set a flag so we generate a warning on screen
                end
                
                if ismember(locX(i_X),locPs) % than we have immobility and death states separately
                    data_out(1,:) = 0; % make all initial things at t=0 zero!
                    % Needed as first row in data for immobility and death
                    % states will be the total amount of animals with which
                    % the experiment started!
                end
                
                means.X{i_X,i_d,i_c} = [t_id data_out];

                % also collect the interpolated model values belonging to the data points!
                loc_c = find(mod_c == mod_c_i_d(i_c)); % where is this concentration in the total model vector from X0mat?
                means.mod{i_X,i_d,i_c} = [t_id interp1(mod_t,XoutX{i_X}(:,loc_c),t_id) interp1(mod_t,XloX{i_X}(:,loc_c),t_id) interp1(mod_t,XhiX{i_X}(:,loc_c),t_id)];
                
                chck_nans = mean(isnan(data_tmp_ic),2); % this is 0 when there are no NaNs, and 1 when there are NaNs in all replicates
                if all(chck_nans == 1 | chck_nans == 0) % if they are all 0 or 1, we can calculate deaths as usual
                    % also collect the predicted deaths from the model survival
                    % and the initial number of animals S0
                    Dmod = means.mod{i_X,i_d,i_c}(:,[2:4]); % extract survival probabilties with CI
                    Dmod = Dmod - [Dmod(2:end,:);zeros(1,3)]; % death probabilities in each interval, from the data
                    means.moddeaths{i_X,i_d,i_c} = [t_id S0.*Dmod(:,1) S0 .*Dmod(:,2) S0.*Dmod(:,3)];
                else % then there is at least one replicate with a NaN(s) and non-NaN(s) at same time point
                    means.moddeaths{i_X,i_d,i_c} = [t_id nan(length(t_id),3)];
                end
                
                means.max(i_X,i_d) = 1; % max(means.max(i_X,i_d),max(min(1,a+b))); % update overal maximum using highest edge of CI
                means.maxt(i_d) = max(means.maxt(i_d),max(t_id)); % update overal maximum using highest edge of CI    
            
            else % then we have length or repro data ...
                
                if transf ~= 1 % then we take mean and se accounting for the transformations (as is done in transfer)
                    lam = 1; % set the transformation indicator to 'none'
                end
                data_out = recalc_data(data_tmp_ic,w_ic,lam);
                plotwts = sum(w_ic,2,'omitnan'); % sum the weights over replicates, take NaNs as zero weight
                % Note: plotwts is only used for plotting outliers
                
                % collect mean and approximate 95% CI of the observation as
                % 2 times the SE
                means.X{i_X,i_d,i_c} = [t_id data_out plotwts];
                
                % also collect the interpolated model values belonging to the data points!
                loc_c = find(mod_c == mod_c_i_d(i_c)); % where is this concentration in the total model vector from X0mat?
                means.mod{i_X,i_d,i_c} = [t_id interp1(mod_t,XoutX{i_X}(:,loc_c),t_id) interp1(mod_t,XloX{i_X}(:,loc_c),t_id) interp1(mod_t,XhiX{i_X}(:,loc_c),t_id)];
                
                means.max(i_X,i_d) = max(means.max(i_X,i_d),max(data_out(:,3))); % update overal maximum using highest edge of CI    
                means.maxt(i_d)    = max(means.maxt(i_d),max(t_id)); % update overal maximum using highest edge of CI    
            end
        end
    end
end

% for i_X = 1:length(locX) % run through trait states (not damage)
%     for i_d = 1:n_D % run through all data sets per state
%         for i_c = 1:length(mod_c_i_d) % run through unique treatments in data set
%             printmat = means.X{i_X,i_d,i_c}; % extract means matrix
%             printmat = printmat(:,1:4); % only keep time, mean and CI
%             if locX(i_X) == 3
%                 printmat(:,2:end) = 1000 * printmat(:,2:end); % lipid times 1000!
%             end
%             
%             if n_D>1 % if there is more than 1 data set
%                 tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
%             else
%                 tit = '';
%             end
%             
%             if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
%                 Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
%                 tit = [tit,' ',Ltmp];
%             else
%                 tit = [tit,' ',leglab1,num2str(mod_c(i_c)),' ',leglab2];
%             end
%             
%             switch locX(i_X)
%                 case 2
%                     tit = [tit,' length (mm)'];
%                 case 3
%                     tit = [tit,' lipid volume (10^-3 mm)'];
%                 case 5
%                     tit = [tit,' oxygen use (umol/d)'];
%             end
%             
%             disp(tit)
%             disp(printmat)
%             disp(' ')
%             
%         end
%     end
% end


%% Plot everything in standard format
% This format is based on the openGUTS plot format, but in a more general
% manner. Now it works for survival data (GUTS), but also sub-lethal traits
% (DEBtox) and binary mixtures.

if flag_meanskip == 1 && plot_repls == 0 % we could not calculate some means, and plot means
    warning('off','backtrace') % no need to display where the warning is generated
    warning('For survival, some means are not plotted because there are NaNs in some replicates.')
    warning('on','backtrace')
end

for i_d = 1:n_D % run through all data sets; make a new figure window for each data set
    
    % mod_c is the complete concentration vector from X0_mat, but we may
    % need less for plotting THIS data set.
    mod_c_i_d = []; % start with empty overall concentration vector (this will collect all c in the data set across all states)
    for i_X = 1:length(locX) % run through all state variables    
        [~,loc_c] = ismember(data.c{i_X,i_d},mod_c); % are the scenarios for this data set in mod_c?
        loc_c(loc_c==0) = []; % remove any zeros (scenario is in data set but not in X0mat)
        mod_c_i_d = cat(2,mod_c_i_d,mod_c(loc_c)); % only that part of mod_c that is used for data set i_d
        % this adds the concentration vector from this state to the complete set!
    end
    mod_c_i_d = unique(mod_c_i_d); % keep unique ones and sort
            
    % assume that the lowest c is the control ... as control will be plotted in all panels
    [~,loc_min] = min(mod_c); % where is the scenario with the minimum overall in model set?
    loc_min = loc_min(1); % just make sure there's only 1 (superfluous)
    
    m      = length(locD) + length(locX); % number of states to plot (number of rows in multiplot)
    n      = length(mod_c_i_d); % number of treatments to plot (columns in multiplot)
    [figh,ft] = make_fig(m,n,2); % create a figure window of correct size
    
    %% Locate maxima for setting y-axis limits for THIS window!
    % Next, I try to find the maximum damage and concentration, ONLY for
    % the treatments plotted in THIS plot window. Otherwise, we get the
    % same scaling across all data sets, which is generally not what we
    % want.
    [~,loc_mod_c]=ismember(mod_c_i_d,mod_c); % find out where the ones we're plotting in this figure are in total mod_c
    for i_D = 1:length(locD) % run through damage states (there may be more than one)
        C_max(i_D) = max(max(XoutC{i_D}(:,loc_mod_c))); % maximum exposure conc. in this plot window
        D_max1 = max(max(XoutD{i_D}(:,loc_mod_c))); % maximum damage level
        D_max2 = max(max(XhiD{i_D}(:,loc_mod_c))); % maximum CI of damage level
        D_max(i_D) = max(D_max1,D_max2); % take maximum of both
    end
    
    % Now do the same for the other state variables.
    % Time vector may differ between data sets, so only use time vector for
    % THIS data set to find maximum for plotting!
    ind_t_max = find(mod_t<=data.maxt(i_d),1,'last');
    
    for i_X = 1:length(locX) % run through all other state variables
        
        if ismember(locX(i_X),locP) % if it is a probability state ...
            X_max(i_X) = 1; % maximum can always be one
            continue % go to next state!
        end

        % First, look for maximum of the model curves
        X_mod_max1 = max(max(XoutX{i_X}(1:ind_t_max,loc_mod_c))); % maximum predicted state level in this plot window
        X_mod_max2 = max(max(XhiX{i_X}(1:ind_t_max,loc_mod_c))); % maximum CI of predicted state level
        X_mod_max  = max([X_mod_max1,X_mod_max2]); % take maximum of both
        if plot_min == 1 % then we also need to look at the controls
            X_mod_max = max([X_mod_max,max(max(XoutX{i_X}(1:ind_t_max,loc_min)))]); % update maximum for control
        end
        
        % Next, look for the maximum of the data (with their CIs) across all data sets
        X_data_max = 0;
        if plot_repls == 1 % then plot individual replicates
            [ind_inset,~]=ismember(data.c{i_X,i_d},mod_c_i_d); % find out where the ones we're plotting in this figure are
            % loc_data_c(loc_data_c==0) = []; % remove any zeros (some states may not have all treatments for this window)
            for i_data_c = 1:length(ind_inset) % run through treatments
                if ind_inset(i_data_c) == 1 % then this treatment will be plotted as it is in mod_c_i_d
                    X_data_max = max([X_data_max , max(max(data.X{i_X,i_d,i_data_c}(:,[2 4])))]);
                end
            end
        else % plot means
            [ind_inset,~]=ismember(means.c{i_X,i_d},mod_c_i_d); % find out where the ones we're plotting in this figure are
            % loc_data_c(loc_data_c==0) = []; % remove any zeros (some states may not have all treatments for this window)
            for i_data_c = 1:length(ind_inset) % run through means of treatments
                if ind_inset(i_data_c) == 1 % then this treatment will be plotted as it is in mod_c_i_d
                    X_data_max = max([X_data_max , max(max(means.X{i_X,i_d,i_data_c}(:,[2 4])))]);
                end
            end
        end
        
        X_max(i_X) = max([X_mod_max,X_data_max]); % keep maximum of model and data
    end
    
    %% Next, move to the plotting part
    for i_cT = 1:n % run through relevant treatments in X0mat
        
        i_c  = find(mod_c == mod_c_i_d(i_cT)); % where is this scenario in the overall set of model curves?   
        row  = -1;
        xlab = glo.xlab; % text on x-axis
        
        % First plot damage and exposure scenario (there are no data for these states!)        
        for i_D = 1:length(locD) % run through damage states (there may be more than one)
            
            row = row + 1; % go to next row of plots for damage
            
            if slowkin == 1 % for slow kinetics
                ylab = 'damage time';
            else 
                ylab = glo.ylab{locD(i_D)}; % text on y axis
            end
            
            mod_plt = [mod_t XoutD{i_D}(:,i_c) XloD{i_D}(:,i_c) XhiD{i_D}(:,i_c)]; % model things to plot: mean, low CI and high CI
            h_pl(row*n+i_cT) = plot_helper(m,n,row*n+i_cT,ft,xlab,ylab,mod_plt,{[]}); % use helper function to plot model curves
            if opt_tktd.plotexp == 1 % if we want to plot the exposure ...
                area(mod_t,XoutC{i_D}(:,i_c),'FaceColor','b','FaceAlpha',0.1); % plot concentration profile as a filled area (option 'FaceAlpha',0.25 removed)
            end
            
            if i_D == 1 && opt_tktd.flip == 0 % only place title in first row
                % create title above panel (but not if we intend to flip later!)
                if n_D>1 % if there is more than 1 data set
                    tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
                else
                    tit = '';
                end
                if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
                    title([tit,Ltmp],ft.name,ft.title)
                else
                    title([tit,leglab1,num2str(mod_c(i_c)),' ',leglab2],ft.name,ft.title)
                end
            end
            
            if lim_data == 1
                xlim([0 data.maxt(i_d)]) % limit x-axis to dataset
            else
                xlim([0 mod_t(end)]) % limit x-axis to model time vector
            end
            if max_on_C == 1
                % ylim([0 max(1e-6,1.05*max(XoutC{i_D}(:)))]) % limit y-axis to exposure level
                ylim([0 max(1e-6,1.05*C_max(i_D))]) % limit y-axis to exposure level
            else
                % ylim([0 max(1e-6,1.05*modmaxD(i_D))]) % limit y-axis to calculated damage
                ylim([0 max(1e-6,1.05*D_max(i_D))]) % limit y-axis to calculated damage level
            end
            
        end
        
        % Plot output states for traits
        for i_X = 1:length(locX) % run through all other state variables
            row = row + 1; % go to next row to plot next state variable
            % collect info for plotting model curve with CI (and plotting
            % the control as well)
            if plot_min == 1
                mod_plt  = [mod_t XoutX{i_X}(:,i_c) XloX{i_X}(:,i_c) XhiX{i_X}(:,i_c) XoutX{i_X}(:,loc_min)];
            else
                mod_plt  = [mod_t XoutX{i_X}(:,i_c) XloX{i_X}(:,i_c) XhiX{i_X}(:,i_c)];
            end
            
            if plot_repls == 1 % then plot individual replicates
                ind_data_c = find(data.c{i_X,i_d} == mod_c(i_c)); % find where this treatment is in the data vector
                % Note: there can be more than one if there are replicated treatments
                clear data_plt % clear it to make sure that no old info remains
                outl_pts = [];
                for i_repl = 1:length(ind_data_c) % run through replicates
                    data_plt{i_repl} = data.X{i_X,i_d,ind_data_c(i_repl)};
                    if ismember(locX(i_X),locP) % then we have quantal data
                        data_plt{i_repl} = data_plt{i_repl}(:,[1:4]); % only keep times, probabilities, and CI lo/hi
                    else % for sub-lethal data ...
                        plotwts = data_plt{i_repl}(:,5); % weights to spot for outliers
                        outl_pts = cat(1,outl_pts,data_plt{i_repl}(plotwts==0,[1 2])); % add outlier points to plot
                    end
                end
                if isempty(ind_data_c) % it is possible that there are no data for this treatment for this state
                    data_plt = {[]}; % just make it empty cell
                end
                ylab = glo.ylab{locX(i_X)}; % label for y-axis
                
                h_pl(row*n+i_cT) = plot_helper(m,n,row*n+i_cT,ft,xlab,ylab,mod_plt,data_plt); % use helper function to plot model and data
                if ~isempty(outl_pts)
                    plot(outl_pts(:,1),outl_pts(:,2),'ko','MarkerSize',10) % plot outliers
                end
                
                if lim_data == 1
                    xlim([0 max(data.maxt(i_d))]) % limit x-axis to dataset
                    ylim([0 X_max(i_X)]); % limit y-axis to data set and model curves for THIS window
                else
                    xlim([0 mod_t(end)]) % limit x-axis to model time vector
                    ylim([0 max(data.max(i_X,i_d),modmaxX(i_X))*1]) % limit y-axis to data set and model curves
                end
                mod_plt = mod_plt(:,2:end); % remove time vector for calculating ylim below
                
                % TEST also plot exposure scenario for internal concentrations
                if ismember(locX(i_X),locC) && opt_tktd.plotexp == 1 % then we have internal concentrations to plot, and want to plot exposure
                    area(mod_t,XoutC{i_D}(:,i_c),'FaceColor','b','FaceAlpha',0.1); % plot concentration profile as a filled area (option 'FaceAlpha',0.25 removed)
                    ylim([0 max(1e-6,max(X_max(i_X),1.05*C_max(i_D)))]) % limit y-axis to exposure level or concentration
                end
                
                % There may be cases where we want to use this function,
                % but without damage state. To make sure that titles are
                % plotted in the first row, we need to repeat this code,
                % which is a bit of a pain.
                if isempty(locD) || opt_tktd.flip == 1
                    if i_X == 1 % only place title in first row
                        % create title above panel
                        if n_D>1 % if there is more than 1 data set
                            tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
                        else
                            tit = '';
                        end
                        if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                            Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
                            title([tit,Ltmp],ft.name,ft.title)
                        else
                            title([tit,leglab1,num2str(mod_c(i_c)),' ',leglab2],ft.name,ft.title)
                        end
                    end
                end
                
            else % we will plot means with CIs
                
                i_c_means = find(means.c{i_X,i_d} == mod_c_i_d(i_cT)); % where is this scenario in the set of mean data? 
                if ~isempty(i_c_means)
                    data_plt{1} = means.X{i_X,i_d,i_c_means}; % take only the relevant means from means.X
                else % it is possible that there are no data for this treatment for this state
                    data_plt{1} = []; % just make it empty
                end
                
                outl_pts = [];
                if ismember(locX(i_X),locP) % then we have quantal data
                    if ~isempty(data_plt{1}) % make sure it is not empty!
                        data_plt{1} = data_plt{1}(:,[1:4]); % only keep times, probabilities with CI
                    end
                elseif ~isempty(data_plt{1}) % ??? ~isempty(outl_pts)
                    plotwts  = data_plt{1}(:,5); % weights to spot for outliers
                    outl_pts = data_plt{1}(plotwts==0,[1 2]); % outlier points to plot
                end
                ylab = glo.ylab{locX(i_X)}; % label for y-axis

                h_pl(row*n+i_cT) = plot_helper(m,n,row*n+i_cT,ft,xlab,ylab,mod_plt,data_plt); % use helper function to plot model and data
                if ~isempty(outl_pts)
                    plot(outl_pts(:,1),outl_pts(:,2),'ko','MarkerSize',10); % plot outlier points
                end
                
                if lim_data == 1
                    xlim([0 max(means.maxt(i_d))]) % limit x-axis to dataset
                    ylim([0 X_max(i_X)]); % limit y-axis to data set and model curves for THIS window
                else
                    xlim([0 mod_t(end)]) % limit x-axis to model time vector
                    ylim([0 max(means.max(i_X,i_d),modmaxX(i_X))*1]) % limit y-axis to data set and model curves
                end
                mod_plt = mod_plt(:,2:end); % remove time vector for calculating ylim below
                
                % There may be cases where we want to use this function,
                % but without damage state. To make sure that titles are
                % plotted in the first row, we need to repeat this code,
                % which is a bit of a pain.
                if isempty(locD)
                    if i_X == 1 % only place title in first row
                        % create title above panel
                        if n_D>1 % if there is more than 1 data set
                            tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
                        else
                            tit = '';
                        end
                        if isfield(glo,'LabelTable') && ismember(mod_c(i_c),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                            Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == mod_c(i_c)}; % look up the label belonging to the j-th scenario
                            title([tit,Ltmp],ft.name,ft.title)
                        else
                            title([tit,leglab1,num2str(mod_c(i_c)),' ',leglab2],ft.name,ft.title)
                        end
                    end
                end
                
            end
        end
    end
    
    % Experimental! Try to swap the first and second row of the plot. This is
    % useful if the plot includes damage AND internal concentrations. Since
    % internal concentrations are part of locX they are plotted AFTER damage,
    % while BEFORE would seem more logical. Code below works, but the titles
    % above the columns now also shift one row down. That is why this option is
    % also used to decide when to plot the titles above. This is quite
    % cumbersome, but other solutions would require more structural revision of
    % the code.
    if opt_tktd.flip == 1
        hFig    = gcf;
        hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
        hAxes   = flipud(hAxes); % reverse order
        for i = 1:n % run through rows of plot
            h_tst = h_pl([n+i i]); % handle for plot on first row and second row, reversed!
            set(hAxes([(i-1)*m+1 (i-1)*m+2]),{'position'},{h_tst.Position}.') % use them to set position in hAxes
            % Note: hAxes has a different order for the plots than h_pl: it
            % goes through the rows rather than through the columns.
        end
    end
    
    % Experimental! Remove y-axis information from first row of plots. This
    % is useful when combining multiple plots.
    if rem_y_info == 1
        figure(figh)
        for yi = 1:m
            axes(h_pl(n*(yi-1)+1))
            set(gca,'YTickLabel',[]);
            set(gca,'YLabel',[])
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
        if n_D == 1
            savenm = ['tktd_plot_',glo.basenm];%
        else
            savenm = ['tktd_plot_set',num2str(i_d),'_',glo.basenm];%
        end
        save_plot(figh,savenm);
    end
    
end


%% Make a predicted-observed plot
% For now, only use means ... at it may be unreadable when plotting
% individual replicates ... Each data set gets a new plot window. I like to
% modify this part such that only states are plotted that actually have
% data ... However, that needs a bit of thinking. Calculated r2 is based on
% the means and includes the resposne at t=0.

if plot_obspred == 0 % unless this option is set to zero ...
    return % then we can just return here
end

n = ceil(sqrt(length(locX)));
m = ceil(length(locX)/n);

for i_d = 1:n_D % run through all data sets
    
    [figh,ft] = make_fig(m,n); % create a figure window of correct size
    
    for i_X = 1:length(locX) % run through all other state variables
        
        h_pl = subplot(m,n,i_X); % make a sub-plot
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        % create axis labels
        xlab = ['obs. ' glo.ylab{locX(i_X)}]; % label
        ylab = ['pred. ' glo.ylab{locX(i_X)}]; % label
        xlabel(xlab,ft.name,ft.label)
        ylabel(ylab,ft.name,ft.label)
        
        % create title above panel
        if n_D>1 % if there is more than 1 data set
            tit = ['Set ',num2str(i_d),', ']; % use set nr. in the title
        else
            tit = '';
        end
        if ismember(locX(i_X),locP) % then we have quantal data
            tit=[tit,'bars: obs. Wilson score, pred. '];
            
            % Plot plus and minus 20% in the plot (roughly equivalent to what
            % is plotted in the EFSA opinion).
            plot([0.2 1],[0 0.8],'k:','LineWidth',1) % plot x% beneath the 1:1 line
            plot([0 0.8],[0.2 1],'k:','LineWidth',1) % plot x% above the 1:1 line
            
        else
            tit=[tit,'bars: obs 2xSE, pred. '];
        end
        
        switch type_conf
            case 0
                tit = [tit,'none'];
            case 1
                tit = [tit,'Bayesian 95% CI'];
            case 2
                tit = [tit,'LikReg 95% CI'];
            case 3
                tit = [tit,'ParSpace 95% CI'];
        end
        
        title(tit,ft.name,ft.title)
        
        mod_c_i_d = means.c{i_X,i_d}; % unique concentration vector in the data set i_d
        
        coll_r_sq = [];
        
        for i_c = 1:length(mod_c_i_d) % run through all modelled treatments for this data set and state
            
            mod_plt = means.mod{i_X,i_d,i_c}(:,[2:4]); % model values with CI
            mod_plt(:,[2 3]) = [mod_plt(:,1)-mod_plt(:,2) mod_plt(:,3)-mod_plt(:,1)]; % recalculate to extent of error bars
            
            data_plt = means.X{i_X,i_d,i_c};
            
            outl_pts = [];
            if ~isempty(data_plt) % make sure it is not empty!
                
                if ~ismember(locX(i_X),locP) % then we have sub-lethal data; any outliers?
                    plotwts  = data_plt(:,5); % weights to spot for outliers
                    outl_pts = [data_plt(plotwts==0,2) mod_plt(plotwts==0,1)]; % outlier points to plot
                end
                
                data_plt = data_plt(:,[2:4]); % only keep state output with CI
                % recalculate CI to extent of the error bars
                data_plt(:,[2 3]) = [data_plt(:,1)-data_plt(:,2) data_plt(:,3)-data_plt(:,1)];
                
            end
            
            errorbar(data_plt(:,1),mod_plt(:,1),mod_plt(:,2),mod_plt(:,3),data_plt(:,2),data_plt(:,3),'ko','MarkerFaceColor','k','LineWidth',1)
            if ~isempty(outl_pts)
                plot(outl_pts(:,1),outl_pts(:,2),'ko','MarkerSize',10); % plot outlier points
            end
            
            max_X = means.max(i_X,i_d);
            plot([0 max_X],[0 max_X],'k:') % plot 1:1 line
            
            coll_r_sq = cat(1,coll_r_sq,[data_plt(:,1) mod_plt(:,1)]); % collect data and model for r-square calculation
            
        end
        
        if ~isempty(mod_c_i_d)
            % R-squared based on MEANS of the data set, and not accounting
            % for transformations! Include response on t=0, which is
            % different from the openGUTS calculations.
            coll_r_sq(isnan(coll_r_sq(:,1)),:) = []; % remove rows where observation is NaN (is that possible?)
            coll_r_sq(isnan(coll_r_sq(:,2)),:) = []; % remove rows where model output is NaN (that is possible when the model time vector is shorter than the data!)
            res     = (coll_r_sq(:,1) - coll_r_sq(:,2)) ; % residuals against model
            res_tot = coll_r_sq(:,1) - mean(coll_r_sq(:,1)) ; % residuals from mean of data
            GoodFit(i_X,i_d) = 1 - (res' * res)/(res_tot' * res_tot); % R-square, NOT adjusted for no. of pars.
        else
            GoodFit(i_X,i_d) = NaN; % this implies that there were no data here
        end
    end
           
    if glo.saveplt > 0 % if we want to save the plot
        if n_D == 1
            savenm = ['tktd_pred_obs_plot_',glo.basenm];%
        else
            savenm = ['tktd_pred_obs_plot_set',num2str(i_d),'_',glo.basenm];%
        end
        save_plot(figh,savenm);
    end
    
end

diary('results.out') % collect output in the diary "results.out"
disp(' ')
disp('Model efficiency (r2) on means, incl. t=0, not accounting for transformations')
disp('================================================================================')
for i_d = 1:n_D % run through all data sets
    if n_D>1
        disp(['Data set ',num2str(i_d)])
    end
    for i_X = 1:length(locX) % run through all other state variables
        if ismember(locX(i_X),locP) % then we have quantal data (there might be more than one state)
            fprintf('  survival    : %#1.4f \n',GoodFit(i_X,i_d))
        elseif ismember(locX(i_X),locL) % then we have length data (there might be more than one state)
            fprintf('  body length : %#1.4f \n',GoodFit(i_X,i_d))
        elseif ismember(locX(i_X),locR) % then we have repro data (there might be more than one state)
            fprintf('  reproduction: %#1.4f \n',GoodFit(i_X,i_d))
        elseif ismember(locX(i_X),locC) % then we have body-residue data (there might be more than one state)
            fprintf('  body residue: %#1.4f \n',GoodFit(i_X,i_d))
        end
    end
    if n_D>1 && i_d<n_D
        disp('================================================================================')
    end
end
disp('================================================================================')
diary off

function h_pl = plot_helper(m,n,i_a,ft,xlab,ylab,mod_plt,data_plt)

% Sub-function to do the plotting and the formatting of the standard plot.

h_pl = subplot(m,n,i_a); % make a sub-plot
hold on
h_ax = gca;
set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
    p = get(h_pl,'position');
    p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
    set(h_pl, 'position', p);
end

if i_a>(n*(m-1)) % only put xlabel on bottom row
    xlabel(xlab,ft.name,ft.label)
else
    set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
end
if i_a==1 || (i_a-1)/n == floor((i_a-1)/n) % only put y labels on first column
    ylabel(ylab,ft.name,ft.label)
else
    set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
end

mod_t = mod_plt(:,1); % read time vector from mod_plt
mod_plt(:,1) = []; % and remove that column

% Little trick to fill the area between the two curves, to
% obtain a coloured confidence interval as a band.
t2  = [mod_t;flipud(mod_t)]; % make a new time vector that is old one, plus flipped one
Xin = [mod_plt(:,2);flipud(mod_plt(:,3))]; % do the same for the plot line, hi and lo
fill(t2,Xin,'g','LineStyle','none','FaceAlpha',1) % and fill this object

plot(mod_t,mod_plt(:,2),'k:','LineWidth',1) % plot dotted line for CI
plot(mod_t,mod_plt(:,3),'k:','LineWidth',1) % plot dotted line for CI

plot(mod_t,mod_plt(:,1),'k-','LineWidth',1) % plot model line for this treatment and this state
% if there is a control line provided, plot it as dotted line
if size(mod_plt,2) > 3 % control will be 4th column
    plot(mod_t,mod_plt(:,4),'k:','LineWidth',1)
end

% if there are data, plot them
if ~isempty(data_plt{1}) % Note that data_plt is always a cell array
    for i_repl = 1:length(data_plt) % run through replicates
        
        if size(data_plt{i_repl},2)>2 % only if there are columns for the CI
            % recalculate CI to extent of the error bars
            errbar = [data_plt{i_repl}(:,2)-data_plt{i_repl}(:,3) data_plt{i_repl}(:,4)-data_plt{i_repl}(:,2)];
            errorbar(data_plt{i_repl}(:,1),data_plt{i_repl}(:,2),errbar(:,1),errbar(:,2),'ko','MarkerFaceColor','k','LineWidth',1)
        else % otherwise, only plot the data point
            plot(data_plt{i_repl}(:,1),data_plt{i_repl}(:,2),'ko','MarkerFaceColor','k')
            error('Does this still happen?') % I think I always fill all columns
        end
        
    end
end
