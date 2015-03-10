function [locs, defs] = circleContainedBrownian(t, param)
% create a random-walk path confined to a circle for timepoints t
%
%  locs = circleContainedBrownian(t, param)
%
%  brownianJitter produces a jitter path based on a brownian motion pattern
%  that resets the location to [0 0] when the distance from [0 0] goes
%  above maxR. The locations are produced for every timepoint t. The
%  parameters provided in the structure param can be used to change
%  jitter-path properties.
%
% param-fields
%   - v         - speed of jitter in x and y in mu/ms -> mm/s       - [1 1]
%   - maxR      - maximum radius before reseting to [0 0]           - 5
%
%  outputs:
%   - locs      - x and y coordinates in 2 columns
%

%% Set defaults
warning('This function is deprecated and not used by current isetbio');
disp('If you really need it, please consider rewrite / check it');
defs.v         = [1 1];
defs.maxR      = 5;

%% Check input
if notDefined('t'), error('Input at least 2 timepoints'); end;
if notDefined('param'), param = struct(); end;

%% Overwrite with user-input
param = setstructfields(defs, param);

%% Extract values;
v    = param.v;         % Speed in x and y axis
maxR = param.maxR;      % Radius of containing circle

%% Calculate jitter path
% Prepare values
dT = [0 t(2 : end) - t(1 : end-1)];     % time between timepoints
nT = length(t);                         % number of timepoints
locs = [0 0];                           % start-location

for iT = 2 : nT
    % Calculate new position
    s  = v * dT(iT);                      % displacement = speed * dT
    dXY = rand(1, 2) .* s - (s / 2);      % add randomness
    locs(iT, :) = locs(iT - 1, :) + dXY;  % new location = old location + displacement
    
    % Check distance from center, and reset if necessary
    if sqrt(sum(locs(iT, :).^2)) > maxR, locs(iT, :) = [0 0]; end;
end