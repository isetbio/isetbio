function current = osCompute(obj, cMosaic, varargin)
% Compute the output response of the L, M and S cone outer segments
%
%   current = osCompute(obj, cMosaic, varargin)
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
% TODO
%  How we are handling the different cone classes is not clear to BW.
%  This needs a scientific review.  For now, to make it run, I am putting
%  in the mean rate for one cone class.  See the code below.
%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% JRG/HJ/BW, ISETBIO TEAM, 2016

%% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osBioPhys'));
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.addParameter('bgR',0,@isnumeric);

% The background absorption rate
p.parse(obj, cMosaic, varargin{:});
bgR = p.Results.bgR;

% This was bgR.  What should it be now?
pRate = cMosaic.absorptions/cMosaic.integrationTime;

%% What should we put in as the bgR in this case?

% JRG and BW need to review the logic here, which may not be right yet.
% How we handle the different cone classes is not clear.
obj.state = osAdaptSteadyState(obj, bgR, [size(pRate, 1) size(pRate, 2)]);

obj.state.timeStep = obj.timeStep;

% How does this handle the separate cone signals?
[current, obj.state]  = osAdaptTemporal(pRate, obj);

% add noise - alert user
if obj.noiseFlag
    disp('Current noise added')
    current = osAddNoise(current);
else
    disp('No current noise added')
end

% In some cases, we run with 1 photoreceptors to set up the LMS filters.
% In that case this is a good curve to plot
%    vcNewGraphWin; plot(squeeze(current));

obj.coneCurrentSignal = current;

end