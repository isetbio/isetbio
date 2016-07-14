function current = osCompute(obj, pRate, coneType, varargin)
% Compute the output response of the L, M and S cone outer segments
%
%   current = osCompute(obj, pRate, coneType, varargin)
% 
% This converts isomerizations (R*) to outer segment current (pA). The
% difference equation model by Rieke is applied here. If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
%
% Inputs: 
%   obj      - the osBioPhys object
%   pRate    - photon absorption rate in R*/sec
%   coneType - cone type matrix, 1 for blank, 2-4 for LMS respectively
% 
% Optional paramters (key-value pairs)
%   'bgR'    - background (initial state) cone isomerization rate
%
% Outputs:
%   current  - outer segment current in pA
%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% JRG/HJ/BW, ISETBIO TEAM, 2016

% check pRate type for backward compatibility
if isstruct(pRate) && isfield(pRate, 'type') && ...
        strcmp(pRate.type, 'sensor')
    warning('The input is a sensor, should update to use coneMosaic.');
    obj.osSet('timestep', sensorGet(pRate, 'time interval'));
    if notDefined('coneType')
        current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
            sensorGet(pRate, 'cone type'));
    else
        current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
            sensorGet(pRate, 'cone type'), coneType, varargin{:});
    end
    % in the old code, we return obj instead of current
    current = obj.osSet('cone current signal', current);
    return
end

% parse inputs
p = inputParser; p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osBioPhys'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);
p.addParameter('bgR', 0, @isscalar);
p.addParameter('append', false, @islogical);

p.parse(obj, pRate, coneType, varargin{:});
bgR = p.Results.bgR;
isAppend = p.Results.append;

% init parameters
if ~isAppend
    p  = osInit;  % default parameters for biophysics model
    obj.state = osAdaptSteadyState(bgR, p, [size(pRate, 1) size(pRate, 2)]);
    obj.state.timeInterval = obj.timeStep;
end

[current, obj.state]  = osAdaptTemporal(pRate, obj.state);

% add noise
if obj.noiseFlag, current = osAddNoise(current); end
if isAppend
    obj.coneCurrentSignal = cat(3, obj.coneCurrentSignal, current);
else
    obj.coneCurrentSignal = current;
end

end