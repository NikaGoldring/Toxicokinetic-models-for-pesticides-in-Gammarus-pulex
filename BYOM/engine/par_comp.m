function par_comp(par,par_tmp,varargin)

% This function compares the _par_tmp_ that it is entered into a certain
% function with the _par_ that is in a saved set.
% 
% Author     : Tjalling Jager 
% Date       : March 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

names_tmp  = fieldnames(par); % extract all field names of par (global)
ind_fittag = ~strcmp(names_tmp,'tag_fitted');
names_tmp  = names_tmp(ind_fittag); % make sure that the fit tag is not in names_tmp
nfields    = length(names_tmp);

pmat_tmp = packunpack(1,par_tmp,0); % transform structure *from input* into a regular matrix
pmat     = packunpack(1,par,0); % transform structure *from saved set* into a regular matrix

warning('off','backtrace')
if size(pmat,1) ~= size(pmat_tmp,1)
    warning('The parameter matrix in the workspace is of a different size than in the saved set! (saved set will be used)')
    warning('off','backtrace'), disp(' ')
    return
end

% Identify parameter(s) to include in the comparison
ind_comp = true(size(pmat,1),1); % no parameters to ignore, unless something was entered
if ~isempty(varargin) 
    set_zero = varargin{1}; % then it is a cell array with the parameter(s) to set to zero
    if ~isempty(set_zero) % we may want to make a parameter zero (esp. background mortality)        
        if ~iscell(set_zero) % for backward compatibility
            set_zero = {set_zero}; % turn it into a cell array with one element
        end
        loc_zero = false(length(names_tmp),1); % create logical index with zeros
        for i = 1:length(set_zero)
            loc_zero = loc_zero == 1 | strcmp(names_tmp,set_zero{i})==1; % add this parameter to loc_zero (logical indexing)
        end
        ind_comp = ~loc_zero; % parameters to use for the comparison
    end
end

prpar = 0; % don't print the parameter set on screen, unless there is reason to
if ~isequal(pmat(ind_comp,:),pmat_tmp(ind_comp,:))
    
    disp(' ')
    if isequal(pmat(ind_comp,1:2),pmat_tmp(ind_comp,1:2)) % ah, it's just the log settings or boundaries
        warning('The log settings and/or boundaries of parameters in the saved set differs from that in the workspace at the moment. This may not hinder the analysis.')
        prpar = 1;
    end
    if ~isequal(pmat(ind_comp,2),pmat_tmp(ind_comp,2))
        warning('The selection of fitted parameters is different in the saved set.')
        prpar = 1;
    end
    if ~isequal(pmat(ind_comp,1),pmat_tmp(ind_comp,1))
        % perhaps add some code here to NOT trigger this warning when the
        % difference is tiny (rounding errors); but take care of parameters
        % that may be zero ...
        
        diff_p = abs(pmat(ind_comp,1)-pmat_tmp(ind_comp,1))./pmat(ind_comp,1); % relative difference
        ind_zero = find(pmat(ind_comp,1)==0 & pmat_tmp(ind_comp,1)~=0); % any zeros that are not the same?
        diff_p(ind_zero) = inf;
        
        if any(diff_p) > 1e-5
            warning('The best-fitting parameters appear to be different from the one in the workspace at the moment. The one in the workspace is used for the best curve, and the set from the saved set is used for the confidence intervals!')
            prpar = 1;
        else
            warning('The best-fitting parameters appear to be different from the one in the workspace at the moment. The one in the workspace is used for the best curve, but the difference is negligible.')
        end
    end
    
    if prpar == 1
        fprintf('Parameter values from saved set \n');
        fprintf('=================================================================================\n');
        
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
end
warning('off','backtrace'), disp(' ')