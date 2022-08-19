function par = startgrid_debtox(par,FEEDB)

% Usage: par = startgrid_debtox(par,FEEDB)
%
% This function prepares a parameter structure for DEBtox analyses, with
% min-max ranges, as required for the parameter-space explorer
% (<calc_parspace>). As input, it needs to know the parameter structure
% <par>. The calculations for openGUTS are taken as basis.
%
% Note: unlike openGUTS, this function will calculate different (tighter)
% bounds for some parameters when <kd> is fixed. For any fixed parameter,
% its range will be set to the best value, which automatically affects
% calculation of the other ranges. For the effect strengths (bs and bb),
% the rules are very much a guess, so be prepared to do some manual tuning
% when the sample runs into a min-max boundary.
% 
% Inputs
% <par>  parameter structure as obtained from the script
% 
% Outputs
% <par>  parameter structure with optimised search ranges in min/max columns
%
% Author     : Tjalling Jager
% Date       : September 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2020, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <startgrid> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under the
% same license as openGUTS (GPLv3). The modifications are only to ensure
% that the code operates in the general BYOM framework.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

global glo glo2 X0mat % make the data set and initial states global variables

%% BLOCK 1. Preliminary things 
% First some preliminary checks. This section looks whether the parameter
% names <kd> and <zb> (or <zs>) are in the parameter vector, which would
% indicate a DEBtox analysis. If they are, their locations in the
% parameter matrix are collected in globals.

names = glo2.names; % names of the model parameters (for capturing slow kinetics)
glo.loc_kd = find(strcmp(names,'kd')==1); 

glo.loc_mi(1) = find(strcmp(names,'zb')==1); % location for energy-budget threshold on which to spot for slow kinetics
glo.loc_mi(2) = find(strcmp(names,'zs')==1); % also collect location for survival threshold for slow kinetics

if isempty(glo.loc_kd) || isempty(glo.loc_mi) % then we have an incorrect parameter definition (for catching slow kinetics)
    error('Was expecting parameters kd and zb (or zs) in the parameter structure ...') 
end

%% BLOCK 2. Find min-max time and exposure concentration across all data sets. 

c_max = 0;   % look for maximum peak concentration across all data sets
c_min = inf; % look for the minimum *non-zero* TWA across all data sets

t = glo.t;
t_max = max(t);

for i_d = 1:size(X0mat,2) % run through all scenarios for calibration
    
    c     = X0mat(1,i_d); % extract scenario identifier (or concentration)
    c_v   = c; % if no exposure profile is specified, simply copy c (it is a constant concentration)
    c_twa = c; % if no exposure profile is specified, c is also the TWA
    if isfield(glo,'int_scen') % if it exists: use it to derive current external conc.
        if ismember(c,glo.int_scen) % is c in the scenario range global?
            c_v = make_scen(-3,c,t); % use make_scen again to derive actual exposure concentration vector
            % the -3 lets make_scen know we are calling for a conc. vector
            c_twa = cumtrapz(t,c_v)/t_max; % CHECK CHECK
            c_twa = c_twa(end);
        end
    end
    c_max = max(c_max,max(c_v));
    c_min = min([c_min,min(c_twa(c_twa>0))]); % look for the minimum *non-zero* TWA across all data sets
    
end

%% BLOCK 3. Find relevant min-max ranges based on the data set.
% Rules of thumb are used to restrict parameter space, based on the data
% set. These rules are specific for DEBtox, and their performance needs to
% be checked carefully.

slowkin = 0; % TEMP to allow forcing thresholds to log-scale

if par.kd(2) == 1 % only if <kd> is fitted
    par.kd = [1 1 0.01 10 0]; % if time is in days, this is a useful range
else % otherwise, make the range equal to the best value
    par.kd = [par.kd(1) 0 par.kd(1) par.kd(1) 1];
end

if par.zb(2) == 1 % only if <zb> is fitted
    par.zb(3) = c_min*(1-exp(-par.kd(3)*(4/24))); % same as for survival ...
    par.zb(4) = 0.99*c_max;
    if FEEDB(1) == 1 && FEEDB(2) == 0 % for this specific combination, damage can be larger than external concentration
        par.zb(4) = 2*c_max; % so increase the threshold
    end
    par.zb(5) = 1; % fit on normal-scale
    if slowkin == 1
        par.zb(5) = 0; % fit on log-scale
    end
    par.zb(1) = round(mean(par.zb([3 4])),2,'significant'); % geomean is in statistics toolbox, but that's not really needed anyway
else % otherwise, make the range equal to the best value
    par.zb = [par.zb(1) 0 par.zb(1) par.zb(1) 1];    
end

if par.zs(2) == 1 % only if <zs> is fitted
    par.zs(3) = c_min*(1-exp(-par.kd(3)*(4/24))); 
    par.zs(4) = 0.99*c_max;
    if FEEDB(1) == 1 && FEEDB(2) == 0 % for this specific combination, damage can be larger than external concentration
        par.zs(4) = 2*c_max; % so increase the threshold
    end
    par.zs(5) = 1; % fit on normal-scale
    if slowkin == 1
        par.zs(5) = 0; % fit on log-scale
    end
    par.zs(1) = round(mean(par.zs([3 4])),2,'significant');
else % otherwise, make the range equal to the best value
    par.zs = [par.zs(1) 0 par.zs(1) par.zs(1) 1];    
end

if par.bs(2) == 1 % only if <bs> is fitted
    % Min is 10% effect at end of test when fast kinetics and threshold
    % zero. Max is that what you gain in damage in half a day, at lowest
    % exposure, at minimum <kd>, is enough to kill you (95% sure) in x hour.
    par.bs(3) = -log(0.9) / (c_max*t_max);
    par.bs(4) = (2^2*0.95) /(par.kd(3)*c_max*exp(-par.kd(3)*t_max*0.5)); % modified from one hour to half a day (the first 2!)
    % par.bw(4) = (24^2*0.95)/(par.kd(3)*c_max*exp(-par.kd(3)*t_max*0.5)); % as in openGUTS
    par.bs(1) = round(mean(par.bs([3 4])),2,'significant');
    par.bs(5) = 0; % fit on log-scale
else % otherwise, make the range equal to the best value
    par.bs = [par.bs(1) 0 par.bs(1) par.bs(1) 1];    
end

if par.bb(2) == 1 % only if <bb> is fitted
    % Min is stress level of 0.x when fast kinetics and threshold zero. Max
    % is that at minimum <kd>, the stress level at the end of the test is
    % still very high. These ranges depend on the moa, and are set with
    % educated guesses ... however, their performance needs to be checked
    % in each case.
    
    if glo.moa(1) == 1 % then we include an effect on assimilation
        par.bb(3) = 0.2 / c_max;
        par.bb(4) = 2 / (c_max * (1-exp(-par.kd(3)*t_max)));
    elseif glo.moa(2) == 1 % then we include an effect on maintenance
        par.bb(3) = 0.2 / c_max;
        par.bb(4) = 10 / (c_max * (1-exp(-par.kd(3)*t_max)));
    elseif glo.moa(3) == 1 % then we include an effect on growth costs
        par.bb(3) = 0.2 / c_max;
        par.bb(4) = 10 / (c_max * (1-exp(-par.kd(3)*t_max)));
    elseif glo.moa(4) == 1 % then we include an effect on repro costs
        par.bb(3) = 0.5 / c_max;
        par.bb(4) = 2000 / (c_max * (1-exp(-par.kd(3)*t_max)));
    elseif length(glo.moa) >= 5 && glo.moa(5) == 1 % then we include an effect on repro hazards
        % for now, use same as maint/growth ...
        par.bb(3) = 0.2 / c_max;
        par.bb(4) = 20 / (c_max * (1-exp(-par.kd(3)*t_max)));
    end
        
    par.bb(1) = round(mean(par.bb([3 4])),2,'significant');
    par.bb(5) = 0; % fit on log-scale
else % otherwise, make the range equal to the best value
    par.bb = [par.bb(1) 0 par.bb(1) par.bb(1) 1];    
end


