function prob = emSteadyState(em, varargin)
% Estimate prior probability distribution of human eye movement position
%
%   prob = emSteadyState(em)
%
% Inputs:
%   em       - eyemovement structure, see emCreate, emGet, emSet
%   varargin - name value pairs for optional arguments, include
%     'n samples': number of samples to be used in prior estimation
%     'sensor': sensor structure to be used in estimation
%
% Outputs:
%   prob  - probability at each cone position
%
% Example:
%   prob = emSteadyState;
%   vcNewGraphWin; imagesc(prob)
%
% Notes:
%   The estimated probability distribution should somehow close to a 2D
%   gaussian distribution. For default parameters, the fitted gaussian is
%   of zero mean and 7.6 cones standard deviation
%
% See also:
%   emCreate, emGet, emSet
%
% HJ/BW, ISETBIO TEAM, 2015

%% Check input
if notDefined('em'), em = emCreate; end

% Init sensor and parameters
sensor = sensorCreate('human');
nSamples = 1e5;

% parse user input
for ii = 1:2:length(varargin)
    switch ieParamFormat(varargin{ii})
        case 'sensor'
            sensor = varargin{ii+1};
        case 'nsamples'
            nSamples = varargin{ii+1};
        otherwise
            error('Unknown parameters')
    end
end

%% Generate samples
sensor = sensorSet(sensor, 'eye movement', em);
sensor = sensorSet(sensor, 'positions', zeros(nSamples, 2));
sensor = emGenSequence(sensor);

%% Summarize eye movement positions
%  get positions
pos = sensorGet(sensor, 'positions');

%  make it symetric about origin
pos = [pos; -pos];

% estimate probability distribution from data
cts = {-40:40, -40:40};
prob = hist3(pos, cts);
prob = prob / sum(prob(:));
end