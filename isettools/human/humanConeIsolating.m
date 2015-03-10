function [coneIsolating, monitor2cones] = humanConeIsolating(dsp)
% Monitor colors (R,G,B) for cone isolation
%
% [coneIsolating,monitor2cones] = humanConeIsolating(displayStructure)
%
% Calculate the three monitor RGB vectors that isolate cone responses for
% the Stockman L,M and S-cones. These RGB values are returned as the
% columns of the 3x3 coneIsolating matrix. The inverse matrix, mapping from
% monitor values into LMS space, is given in monitor2cones.
%
% The units used in these calculations are designed around the normalized
% RGB values from the monitor space.  Hence, the LMS values are not in
% physical coordinates (e.g., absorptions).  
%
% To calculate this use the display spectral power distribution of the
% monitor, which is energy units of watts/sr/...
%
% Example:
%   dsp = displayCreate('LCD-Apple');
%   signalDirs = humanConeIsolating(dsp)
%   coneIsolatingSPD = displayGet(dsp,'spd')*signalDirs; 
%   vcNewGraphWin; plot(displayGet(dsp,'wave'),coneIsolatingSPD);
%
% Copyright ImagEval Consultants, LLC, 2005.

warning('Deprecated, use displayGet(d, ''rgb2lms'') instead');
if notDefined('dsp'), dsp = displayCreate;  end

wave = displayGet(dsp,'wave');
cones = ieReadSpectra('stockman',wave);  % Energy

spd = displayGet(dsp,'spd');             % Energy

% cones*spd*(r,g,b)
monitor2cones = cones'*spd;

% First column of coneIsolating, col1, maps into cone isolating direction
%
%    monitor2cones * col1 = (1,0,0)
%
% and so forth 
if size(monitor2cones,1) == size(monitor2cones, 2)
    coneIsolating = inv(monitor2cones);
else
    warning('monitor to cones not invertible, returning empty');
    coneIsolating = [];
end