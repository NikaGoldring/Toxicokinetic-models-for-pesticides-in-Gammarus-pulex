%% BYOM function pathdefine.m (setting the path to the engine)
%
%  Syntax: pathdefine
%
% A tiny piece of 'intelligent' code to add the required directories to the
% Matlab Path. Works both for UNIX and PC, as long as the byom script from
% which this function is called is in a subdirectory of the BYOM directory.
% The addition to the path is temporary, and forgotten when Matlab is
% restarted.
%
% * Author: Tjalling Jager 
% * Date: May 2020
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function pathdefine

rem    = pwd; % remember previous location
curdir = pwd; % current directory as string

%% Search and add
% This section goes up the path from the present directory until it has
% located the BYOM root directory. From there, it adds the engine directory
% to the path and returns to the original directorty. This works both on PC
% and Unix (and probably also on Linux)

if isempty(findstr(['BYOM',filesep,'engine'],path)) % BYOM/engine is not already in the path ... add it!
    while strcmp(curdir(end-3:end),'BYOM')==0 % as long as we are not yet in the BYOM folder ...
        cd ..; % go one directory up
        curdir = pwd; % current directory as string
        if length(curdir) < 4 % if dir name is very small, we are probably in the root
            eval(['cd(',char(39),rem,char(39),')']); % go back to the original location
            error('Pathdefine.m cannot find the directory BYOM. Make sure your script is in a sub-directory of BYOM.')
        end
    end
    addpath([curdir,filesep,'engine'],'-begin'); % add engine to the path
    eval(['cd(',char(39),rem,char(39),')']); % go back to the original location
    % This has to be done in this rather complicated manner as windows file locations may include spaces, 
    % which cd interprets in the wrong way!
end

% For now do the same thing again for the parspace directory (this can be
% done smarter, but I would like to avoid problems with older path_define
% versions and BYOM versions without parspace in the engine)
if isempty(findstr(['BYOM',filesep,'parspace'],path)) % BYOM/parspace is not already in the path ... add it!
    while strcmp(curdir(end-3:end),'BYOM')==0 % as long as we are not yet in the BYOM folder ...
        cd ..; % go one directory up
        curdir = pwd; % current directory as string
        if length(curdir) < 4 % if dir name is very small, we are probably in the root
            eval(['cd(',char(39),rem,char(39),')']); % go back to the original location
            error('Pathdefine.m cannot find the directory BYOM. Make sure your script is in a sub-directory of BYOM.')
        end
    end
    if exist([curdir,filesep,'engine',filesep,'parspace']) == 7 % the directory parspace exists
        addpath([curdir,filesep,'engine',filesep,'parspace'],'-begin'); % add parspace to the path
    end
    eval(['cd(',char(39),rem,char(39),')']); % go back to the original location
    % This has to be done in this rather complicated manner as windows file locations may include spaces, 
    % which cd interprets in the wrong way!
end