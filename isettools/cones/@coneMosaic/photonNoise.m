function [noisyImage, theNoise, seed] = photonNoise(absorptions,varargin)
%PHOTONNOISE  Add photon noise to the absorptions 
%   [noisyImage, theNoise, seed] = photonNoise(absorptions,varargin)
%
%   Photon noise is Poisson.  The Poisson variance is equal to the mean. We trap
%   the cases when the value is small (less than 25) and use
%   real Poisson random value, which is slower to compute. At
%   levels above lambda = 25, we approximate the noise as Gaussian.
%
%   The decision to add noise is set in the coneMosaic.noiseFlag.  We are
%   only here if noiseFlag is 'random' or 'frozen'
%
%   Inputs:  
%   absorptions - typically coneMosaic.absorptions
%
%   Outputs:
%   noiseImage - absorptions plus noise
%   theNoise   - the noise that was added to the absorptions
%   seed       - The rng(seed) we used
%
%   Optional parameter name/value pairs chosen from the following:
%
%   'noiseFlag'       Specify frozen noise ('frozen') or not ('random', default)
%   'seed'            When using frozen noise, you can set a seed (default 1)

% HJ/BW ISETBIO Team 2016

%%
p = inputParser;
p.addRequired('absorptions',@isnumeric); 

vFunc = @(x)(ismember(x,{'random','frozen'}));
p.addParameter('noiseFlag','random',vFunc);
p.addParameter('seed',1,@isnumeric);      % Seed for frozen noise

p.parse(absorptions,varargin{:});

% absorptions  = p.Results.absorptions;
seed         = p.Results.seed;
noiseFlag    = p.Results.noiseFlag;

%% Set up RNG depending on noiseFlag

switch noiseFlag
    case 'frozen'
        rng(seed);
    case 'random'
        rng('shuffle');
end

%% This is std * N(0,1)
theNoise = sqrt(absorptions) .* randn(size(absorptions));

% We add the mean electron and noise electrons together.
noisyImage = round(absorptions + theNoise);

% Now, we find the small mean values and create a true Poisson sample. The
% Poisson algorithm is slow for big numbers, but it is fast enough for
% small numbers. We can't rely on the Stats toolbox being present, so we
% use this Poisson sampler from Knuth.

% Apply the Poisson when the mean is less than this
poissonCriterion = 25;
idx = find(absorptions(:) < poissonCriterion);
v = absorptions(absorptions < poissonCriterion);
if ~isempty(v)
    % If noiseFlag is 'random', this routine ignores the seed.
    vn = iePoisson(v,'noiseFlag',noiseFlag,'seed',seed);  % Poisson samples
    noisyImage(idx) = vn;
    if nargout > 1
        % Saves time
        theNoise(idx) = vn - absorptions(idx);
    end
end

% Find the highly unusual case in which the sum of the mean and noise are
% less than zero. This only happens if the mean is 25 or greater and the
% noise is less than -25.  Very unusual, but it can happen given how many
% times we run this random number generate.  So this little clipping is an
% imperfection but it isn't the worst thing we do.
idx = (noisyImage < 0);
noisyImage(idx)  = 0;
if nargout > 1
    % Saves time
    theNoise(idx)    = noisyImage(idx) - absorptions(idx);
end

end