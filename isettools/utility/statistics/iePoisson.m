function [val, seed] = iePoisson(lambda, varargin)
% Create a matrix of Poisson samples using rate parameters in lambda
%
%    [val, seed] = iePoisson(lambda,'nSamp',N,'noiseFlag',...,'seed',...)
%
% The lambda (rate parameter) can be a scalar and you can request multiple
% samples (nSamp) or it can be a matrix of rate parameters.
%
% The noiseflag {'random' or 'frozen'} determines whether to generate a new
% noise seed or to set the noise seed to a fixed value (default is 1).
%
% We first try to call the mex function iePoissonrnd.
% If this does not exist, we try to use the statistics toolbox
% If this fails, we use an algorithm  from Knuth.
%
% Examples (see Run Configuration in Debug):
% % Matrix form
%   nSamp = 128; lambda = round(rand(nSamp,nSamp)*10);
%   tic, val = iePoisson(lambda); toc
%   vcNewGraphWin;
%   subplot(2,1,1), imagesc(lambda); colormap(gray);colorbar; axis image
%   subplot(2,1,2), imagesc(val); colormap(gray); colorbar; axis image
%
% % Multiple samples form
%   lambda = 4; nSamp = 1000; val = iePoisson(lambda,'nSamp',nSamp);
%   vcNewGraphWin; hist(val,50)
%
% % Frozen noise
%     lambda = 4; nSamp = 1; 
%     val1 = iePoisson(lambda,'nSamp',nSamp,'noiseFlag','frozen');
%     val2 = iePoisson(lambda,'nSamp',nSamp,'noiseFlag','frozen');
%     val1, val2
%
% % Return seed
%     [val,seed] = iePoisson(lambda,'nSamp',nSamp,'noiseFlag','random');
%
% Reference: Knuth routine - found on a web-page reference and also at
%   http://en.wikipedia.org/wiki/Poisson_distribution
%   http://www.columbia.edu/~mh2078/MCS04/MCS_generate_rv.pdf
%
% Copyright ImagEval, LLC, 2010
%
% 6/3/15  xd  iePoissrnd now uses a randomly generated seed
% 6/4/15  xd  added flag to determine if noise should be frozen
% 1/24/17 npc Now checking first for Statistics Toolbox, and using 
%             poissrnd if it exists. If not it uses the local equivalent,
%             iePoissrnd

%% Parse parameters
p = inputParser;

% Required
p.addRequired('lambda',@isnumeric);

% Key/value
p.addParameter('nSamp',1,@isnumeric);
vFunc = @(x)(ismember(x,{'random','frozen'}));
p.addParameter('noiseFlag','random',vFunc);
p.addParameter('seed',1,@isnumeric);

p.parse(lambda,varargin{:});

nSamp     = p.Results.nSamp;
noiseFlag = p.Results.noiseFlag;
seed      = p.Results.seed;

if checkToolbox('Statistics and Machine Learning Toolbox')
    switch noiseFlag
        case 'frozen'
            rng(seed);
        case 'random'
            seed = rng;
    end

    % Matlab toolbox version is present. Use it.
    if isscalar(lambda)
        val = poissrnd(lambda, nSamp);
    else
        val = poissrnd(lambda);
    end
       
    return;
end

%% Check if we have MEX function
if (exist('iePoissrnd','file')==3)
    switch noiseFlag
        case 'frozen'
            val = iePoissrnd(lambda, nSamp, seed);
        case 'random'
            seed = rand * 12345701;
            val = iePoissrnd(lambda, nSamp, rand * 12345701);
    end
    return;
end



%% No toolbox, no mex file, making do with home grown Knuth code
switch noiseFlag
    case 'frozen'
        rng(seed);
    case 'random'
        seed = rng;
end
    
% Use the local method. This is slow for large lambda, but a desperate
% move.
warning('Using slow Poission random variable generation');
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

end