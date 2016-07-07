function val = iePoisson(lambda, nSamp, newNoise)
% Create a matrix of Poisson samples using rate parameters in lambda
%
%   val = iePoisson(lambda,nSamp,[newNoise])
%
% The rate parameter can be a scalar, requesting multiple samples, or it
% can be a matrix of rate parameters.
%
% The newNoise flag determines whether to generate a new noise seed or to
% set the noise seed to a fixed value (1).  If newNoise is false, we are
% using frozen noise with a seed = rng(1);
%
% This algorithm is from Knuth.
%
% This script replaces the Matlab poissonrnd function because that function
% is only in the stats toolbox.
%
% The function is used in ISET when we find the pixels with a mean of less
% than, say 10, in which case the normal approximation doesn't work well.
% In that case we over-write the Gaussian shot noise at those locations
% with  Poisson random values. We could use this for all of the values if
% it turns out this routine runs fast enough.
%
% See also: noiseShot
%
% Examples (see Run Configuration in Debug):
% % Matrix form
%   nSamp = 128;
%   lambda = round(rand(nSamp,nSamp)*10);
%   tic, val = iePoisson(lambda); toc
%   figure(1); clf
%   subplot(2,1,1), imagesc(lambda); colormap(gray);colorbar; axis image
%   subplot(2,1,2), imagesc(val); colormap(gray); colorbar; axis image
%
% % Multiple samples form
%   lambda = 4; nSamp = 1000;
%   val = iePoisson(lambda,nSamp);
%   figure(1); clf, hist(val,50)
%
% Reference: Knuth routine - found on a web-page reference and also at
%   http://en.wikipedia.org/wiki/Poisson_distribution
%   http://www.columbia.edu/~mh2078/MCS04/MCS_generate_rv.pdf
%
% Copyright ImagEval, LLC, 2010
%
% 6/3/15  xd  iePoissrnd now uses a randomly generated seed
% 6/4/15  xd  added flag to determine if noise should be frozen

if notDefined('lambda'), error('rate parameter lambda required'); end
if notDefined('nSamp'), nSamp = 1; end
if notDefined('newNoise'), newNoise = true; end

% Check if we have MEX function
if (exist('iePoissrnd','file')==3)
    if newNoise  % Create a new seed and use it
        val = iePoissrnd(lambda, nSamp, rand * 12345701);
    else        % No new noise seed
        val = iePoissrnd(lambda, nSamp);
    end
    return;
end

% Check for stats toolbox
if checkToolbox('Statistics Toolbox')
    
    % Set the RNG to this seed, reliably
    if ~newNoise, p = rng; rng(1); end
    
    % Matlab toolbox version is present. Use it.
    if isscalar(lambda)
        val = poissrnd(lambda, nSamp);
    else
        val = poissrnd(lambda);
    end
    
    % Return the RNG to the current state
    if ~newNoise, rng(p); end
    
    return 
end


% Frozen noise case.  Use sd = 1;
if ~newNoise, p = rng; rng(1); end

% Use the local ISET methods
% Simple implementation, this is slow for large lambda
% Not recommanded, should try to use mex file first
warning('Using slow poission random variable generation');
if isscalar(lambda)
    % Scalar version of the routine
    % Probably we want multiple samples for a single lambda
    val = zeros(1, nSamp);
    for nn = 1 : nSamp
        kk = 0;
        p = 0;
        while p < lambda
            kk = kk + 1;
            p = p - log(rand);
        end
        val(nn) = kk - 1;
    end
    % figure(1); hist(val,50)
else
    % A matrix or vector of lambdas and we return samples for each
    val = zeros(size(lambda));
    for ii = 1 : numel(lambda)
        kk = 0;
        p = 0;
        while p < lambda(ii)
            kk = kk + 1;
            p = p - log(rand);
        end
        val(ii) = kk - 1;
    end
end

% Put the system back
if ~newNoise, rng(p); end
end