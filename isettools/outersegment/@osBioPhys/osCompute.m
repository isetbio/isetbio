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

% parse inputs
p = inputParser; p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osBioPhys'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);
p.addParameter('bgR', 0, @isscalar);

p.parse(obj, pRate, coneType, varargin{:});
bgR = p.Results.bgR;

% init parameters
p  = osInit;  % default parameters for biophysics model

initialState = osAdaptSteadyState(bgR, p, [size(pRate, 1) size(pRate, 2)]);
initialState.timeInterval = obj.timeStep;

current  = osAdaptTemporal(pRate, initialState);

% add noise
if obj.noiseFlag, current = osAddNoise(current); end
obj.osSet('coneCurrentSignal', current);

end