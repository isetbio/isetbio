function clx
% Script to clear many things
%
% Syntax:
%   clx;
%
% Description:
%    Clear many things, including: mex, all, fun
%       and
%    Close many things, including: all, hidden
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    06/19/03  ARW  Created
%    12/12/17  JNM  Formatted
%    01/19/18  jnm  Formatting update to match Wiki.

clear mex
clear all

% Important to close windows prior to clearing variables
% (sometimes the closerequestfunction on a window
% is set to something funny -- want to force it
% to close:)
figs = get(0,'Children');
set(figs,'CloseRequestFcn','closereq');
close all
if exist('hidden', 'var')
    close hidden
end

clear all
clear fun
