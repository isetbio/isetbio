%% ieInit
%
% This script 
%
%   * Closes down current instance of ISETBIO
%   * Closes figures
%   * Starts a fresh version of ISETBIO
%
%   * User variables are NOT cleared
%
% BW, ISETBIO Team, Copyright 2015

%% Check if there is a running instance, closes it, and clears workspace
ieMainClose;   % Closes open ISETBIO windows
close all;     % Closes other windows (e.g., vcNewGraphWin)]
clear global   % Clear out pesky globals

%% Initialize variables

global vcSESSION; %#ok<NUSED>

thisVersion = 4.0; 
ieSessionSet('version',thisVersion);

%% Default session file name.

% There was a thought back in the day that people would want to preserve
% their computed data, and we would store the computations in a session
% file.  This is the code that sets up the name for the current sessions.
%
% In practice, it seems that people never do this, they just recompute.
%
% We check for a session file named iset-dateTime
%
%   * If one exists, we load it.  
%   * If several exist, we load the latest one
%   * The user can load a sessison file with a different name from the Main
%     Window.
%

ieInitSession;           % Initializes the vcSESSION database variables
ieSessionSet('dir',pwd);
sessionFileName = 'isetSession.mat';
ieSessionSet('name',sessionFileName);

%% Clean up, go home
clear sessionFileName
clear thisVersion

%% END
