function [rnd,par,par_sel] = load_rnd(opt_conf)

% This function handles the loading of the sample from parameter space as
% saved by <calc_slice>, <calc_likreg> and <parspace>. A sub-sample
% can/will be used, which is handled here as well. The output <par_sel> may
% not be used ... The options structure <opt_conf> (see <prelim_checks>) is
% used here as well.
%
% Author     : Tjalling Jager 
% Date       : March 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

filenm = glo.basenm;

type_conf = opt_conf.type;     % use values from slice sampler (1), likelihood region(2) to make intervals
lim_set   = opt_conf.lim_set;  % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
n_lim     = opt_conf.n_lim;    % size of limited set (likelihood-region only)
crit_add  = opt_conf.crit_add; % small addition to chi2 criterion for inner rim (to make sure coverage is adequete)
set_zero  = opt_conf.set_zero; % parameter name(s) to set to zero (usually the background hazard)

% for the sample of the likelihood-based joint confidence region, we make a
% small addition to chi2 criterion for the inner rim (to make sure coverage
% is adequete, since we use a discrete, limited, sample)
chicrit   = 3.8415+crit_add; % critical value for 95% confidence at df=1, with an addition
chicrit2  = 3.8415-2*crit_add; % with a subtraction ...

disp(' ')
switch type_conf
    case 1 % sample from posterior distribution, take 95% of MODEL CURVES
        if exist([filenm,'_MC.mat'],'file') ~= 2
            disp('There is no MCMC sample saved, so for intervals run calc_slice first!')
            rnd = -1;
        else
            load([filenm,'_MC'],'rnd','par','par_sel') % load the random sample from the last MCMC run
            % this loads the random sample in rnd, and the parameter matrix par and the
            % selection matrix par_sel
            if lim_set == 2
                disp('  There is no possibility to use an outer hull for the Bayesian analysis.')
                disp(['  Full sample of ',num2str(size(rnd,1)),' sets used.'])
            elseif lim_set == 1 && n_lim < size(rnd,1) % but we can take a sub-sample! (cannot take more than there are)
                rnd = rnd(randperm(size(rnd,1),n_lim),:); % take n_lim random sets from outer hull
                disp(['  Sub-sample of ',num2str(size(rnd,1)),' random sets used.'])
            else
                disp(['  Full sample of ',num2str(size(rnd,1)),' sets used.'])
            end
            disp('Calculating CIs on model predictions, using posterior from slice sampler ... please be patient.')
            
            % rnd is the full sample, with a minlog-likelihood value in the last column
            rnd = rnd(:,1:end-1); % remove the likelihood column (probably not filled anyway)
        end
        
    case 2 % LH sample from joint conf. region, take min and max of curves
        if exist([filenm,'_LR.mat'],'file') ~= 2
            disp('There is no confidence set saved, so for intervals run calc_likregion first!')
            rnd = -1;
        else
            load([filenm,'_LR'],'rnd','par','par_sel') % load the random sample from the last calc_likregion run
            % this loads the sets from the conf region in rnd, and the
            % parameter matrix par and the selection matrix par_sel. Note
            % that rnd is already sorted based on likelihood.
            % 
            % last column in rnd is 2 times loglikratio
            
            rnd = rnd(rnd(:,end)<chicrit,:); % keep the sets in inner rim (df=1) 
            % note that we have added a little bit to the chi2 criterion as we have a
            % limited, discrete sample
            if crit_add > 0
                disp(['Amount of ',num2str(crit_add),' added to the chi-square criterion for inner rim'])
                disp('Slightly more is taken to provide a generally conservative estimate of the CIs.')
            end
            switch lim_set
                case 0 % use full set within inner rim
                    disp(['  Full sample of ',num2str(size(rnd,1)),' sets used.'])
                case 1 % use random sub-sample from out hull (fast)
                    rnd = rnd(rnd(:,end)>chicrit2,:); % keep the sets in outer part of inner rim (df=1)
                    if n_lim > size(rnd,1) % cannot take more than there are
                        disp(['  Outer hull of ',num2str(size(rnd,1)),' sets used (not enough points for a more limited set).'])
                    else
                        rnd = rnd(randperm(size(rnd,1),n_lim),:); % take n_lim random sets from outer hull
                        disp(['  Limited set of ',num2str(size(rnd,1)),' random sets from outer hull used.'])
                    end
                case 2 % use complete outer hull ...
                    rnd = rnd(rnd(:,end)>chicrit2,:); % keep the sets in outer part of inner rim (df=1)
                    disp(['  Outer hull of ',num2str(size(rnd,1)),' sets used.'])
            end
            rnd = rnd(:,1:end-1); % and remove last column
            disp('Calculating CIs on model predictions, using sample from joint likelihood region ... please be patient.')
        end
        
    case 3 % for parameter space explorer
        if exist([filenm,'_PS.mat'],'file') ~= 2
            disp('There is no confidence set saved, so for intervals run optimisation with parameter-space option first!')
            rnd = -1;
        else
            load([filenm,'_PS'],'pmat','coll_all') % load the random sample from the last parameter space run
            % this loads the sets from the conf region in coll_all, and the
            % parameter matrix pmat. Note that coll_all is already sorted based on likelihood.
            %
            % last column in coll_all is log-likelihood itself
            
            par_sel = pmat(:,2); % this is the selection vector
            % put parameters that need to be fitted on log scale back on normal scale
            pmat(pmat(:,5)==0,1) = 10.^(pmat(pmat(:,5)==0,1));
            par = packunpack(2,[],pmat); % pack pmat into structure par
            
            mll = coll_all(1,end); % extract maximum likelihood value from saved set (<coll_all> is sorted, so it is in the first row)
            rnd = coll_all; % rename coll_all to rnd
            clear coll_all; % clear coll_all to save memory
            
            rnd(:,end) = 2 * (rnd(:,end)-mll); % translate into 2 times loglik ratio (as for other samples)            
            
            rnd = rnd(rnd(:,end)<chicrit,:); % keep the sets in inner rim (df=1) 
            % note that we have added a little bit to the chi2 criterion as we have a
            % limited, discrete sample
            if crit_add > 0
                disp(['Amount of ',num2str(crit_add),' added to the chi-square criterion for inner rim'])
                disp('Slightly more is taken to provide a generally conservative estimate of the CIs.')
            end
            switch lim_set
                case 0 % use full set within inner rim
                    disp(['  Full sample of ',num2str(size(rnd,1)),' sets used.'])
                case 1 % use random sub-sample from out hull (fast)
                    rnd = rnd(rnd(:,end)>chicrit2,:); % keep the sets in outer part of inner rim (df=1)
                    if n_lim > size(rnd,1) % cannot take more than there are
                        disp(['  Outer hull of ',num2str(size(rnd,1)),' sets used (not enough points for a more limited set).'])
                    else
                        rnd = rnd(randperm(size(rnd,1),n_lim),:); % take n_lim random sets from outer hull
                        disp(['  Limited set of ',num2str(size(rnd,1)),' random sets from outer hull used.'])
                    end
                case 2 % use complete outer hull ...
                    rnd = rnd(rnd(:,end)>chicrit2,:); % keep the sets in outer part of inner rim (df=1)
                    disp(['  Outer hull of ',num2str(size(rnd,1)),' sets used.'])
            end
            rnd = rnd(:,1:end-1); % and remove last column
            disp('Calculating CIs on model predictions, using sample from parspace explorer ... please be patient.')
        end
        
    otherwise        
        error('Other optimisation methods are not implemented anymore/yet.')
                
end

% Since I am now using this function in several cases to extract only <par>
% from the saved file, it is handy to also apply the <set_zero> option.
% That way we can easily make predictions for long exposure profiles
% without needing to set <hb> to zero manually. In several functions, this
% setting to zero is also done, such as for the calculations of LCx and LPx
% values. So, that means some duplication.
if ~isempty(set_zero) % we may want to make a parameter zero (esp. background mortality)
    % Do not make entire random sample zero as it may hurt the sensitivity analysis
    
    % allow parameters to be set to zero, such as background hazard or initial concentrations
    if ~iscell(set_zero) % for backward compatibility
        set_zero = {set_zero}; % turn it into a cell array with a single element
    end
    for i = 1:length(set_zero)
        eval(['par.',set_zero{i},'(1) = 0;']); % set parameter to zero in par
    end
    
end
