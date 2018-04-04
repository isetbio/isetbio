function current = osCompute(obj, cMosaic, varargin)
% Compute the photocurrent response of the L, M and S cones
%
% Syntax:
%   current = osCompute(obj, cMosaic, varargin)
%
% Description:
%    This converts isomerizations (R*) to outer segment current (pA). The
%    difference equation model by Rieke is applied here. 
%
%    If the noiseFlag property of the osLinear object is set to true, this
%    method adds photocurrent noise to the output signal. See osAddNoise().
%
% Inputs: 
%    obj      - the osBioPhys object
%    cMosaic  - The parent of the outersegment object
%
% Outputs:
%   current   - outer segment current in pA
% 
% Optional key/value pairs:
%    seed     - Seed to use when obj.noiseFlag is 'frozen' (default = 1)
%    'bgR'    - background (initial state) cone isomerization rate
%
% References:
%    http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%    https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%

% History:
%    xx/xx/16  JRG/HJ/BW  ISETBIO TEAM, 2016
%    02/14/18  jnm        Formatting

%% parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osBioPhys'));
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.addParameter('bgR', 0, @isnumeric);
p.addParameter('seed', 1, @isnumeric);

% The background absorption rate
p.parse(obj, cMosaic, varargin{:});

seed = p.Results.seed;

% This is the background isomerization rate mean for each cone class
% (R*/sec). It could be calculated here using coneMeanIsomerizations. Not
% sure why we pass it in.
bgR = mean(p.Results.bgR);

% R*/sec over time (x, y, t) for each one. 
pRate = cMosaic.absorptions / cMosaic.integrationTime;

%% What should we put in as the bgR in this case?
% [TODO: JRG and BW need to review the logic here, which may not be right
% yet. How we handle the different cone classes is not clear.]
obj.state = osAdaptSteadyState(obj, bgR, [size(pRate, 1) size(pRate, 2)]);

obj.state.timeStep = obj.timeStep;

% How does this handle the separate cone signals?
[current, obj.state] = osAdaptTemporal(pRate, obj);

% The outer segment noise flag
if (~strcmp(obj.noiseFlag, 'none'))  
    switch obj.noiseFlag
        case 'frozen'
            rng(seed);
        case 'random'
            rng('shuffle');
    end
    current = osAddNoise(current);
    fprintf('Added osBiophys noise (''%s'')', obj.noiseFlag);
end

end