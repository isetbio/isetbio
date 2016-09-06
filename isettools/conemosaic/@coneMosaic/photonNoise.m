
function [noisyImage, theNoise] = photonNoise(absorptions,varargin)
% Photon noise at the absorptions is Poisson.
% The Poisson variance is equal to the mean
% For moderate mean levels, Poisson is very close to Normal
%
% We multiply each point in the absorption data by the square
% root of its mean value to create the noise standard
% deviation. For most cases the Normal approximation is
% adequate. We trap (below) the cases when the value is small
% (less than 25) and replace it with the real Poisson random
% value, which is slower to compute.

p = inputParser;
p.addRequired('absorptions',@isnumeric);

% Frozen or new random noise
p.addParameter('newNoise',true,@islogical);

p.parse(absorptions,varargin{:});
absorptions  = p.Results.absorptions;
newNoise = p.Results.newNoise;

% This is std * N(0,1)
theNoise = sqrt(absorptions) .* randn(size(absorptions));

% We add the mean electron and noise electrons together.
noisyImage = round(absorptions + theNoise);

% Now, we find the small mean values and create a true Poisson
% sample. The Poisson algorithm is slow for big numbers, but it
% is fast enough for small numbers. We can't rely on the Stats
% toolbox being present, so we use this Poisson sampler from
% Knuth.

% Apply the Poisson when the mean is less than this
poissonCriterion = 25;
idx = find(absorptions(:) < poissonCriterion);
v = absorptions(absorptions < poissonCriterion);
if ~isempty(v)
    vn = iePoisson(v, 1, newNoise);  % Poisson samples
    theNoise(idx) = vn - absorptions(idx);
    noisyImage(idx) = vn;
end

% Find the highly unusual case in which the sum of the mean and
% noise are less than zero. This only happens if the mean is 25
% or greater and the noise is less than -25.  Very unusual, but
% it can happen given how many times we run this random number
% generate.  So this little clipping is an imperfection but it
% isn't the worst thing we do.
idx = (noisyImage < 0);
noisyImage(idx)  = 0;
theNoise(idx)    = noisyImage(idx) - absorptions(idx);
end