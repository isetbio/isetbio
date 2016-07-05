function ir = irCompute(ir, input, varargin)
% Computes the rgc mosaic responses to an input
%
%    ir = irCompute(ir, input, varargin)
%
% The input can be of several types, for historical reasons.
%
% Inputs:
%   ir: inner retina object,
%   input - There are various type of possible input objects.
%
%      'osDisplayRGB'
%      'osIdentity'
%      'osLinear'
%      'osBioPhys'
%      'bipolar'
%
% Computes the responses for all the mosaics attached to the inner retina
% object.
%
% For each mosaic, the linear response is computed by spatial convolution
% of the center and surround RFs. Then, the temporal responses for the
% center and surround is computed, separably. This stage of the computation
% is stored in responseLinear.
%
% The spikes are computed irComputeSpikes routine.  Each rgc type has its
% own compute method.  
%
% Outputs:
%  ir: the inner retina object with responses attached to each mosaic
%
% Example:
%   ir.compute(identityOS);
%   irCompute(ir,identityOS);
%
% See also: rgcGLM, irComputeSpikes,
%
% JRG (c) isetbio team

%% PROGRAMMING TODO
% @JRG
% We plan to move the linear computation to a separate routine,
% irComputeContinuous (or Linear) and carry out the first part of this code
% in that routine, rather than as the switch statement that is here now
%

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));

vFunc = @(x) ~isempty(validatestring(class(x),...
    {'osDisplayRGB','osIdentity','osLinear','osBioPhys','bipolar'}));
p.addRequired('input',vFunc);

p.parse(ir,input,varargin{:});

%% Linear stage of the computation
ir = irComputeLinearSTSeparable(ir, input);

%% Compute spikes for each trial
switch class(ir.mosaic{1})
    case {'rgcLinear'};
        % No linear response implemented yet.
        error('Not yet implemented');
    case{'rgcPhys'}
        % See sublcass routine @irPhys/irCompute.m
    otherwise
        % Runs for rgcLNP, rgcGLM
        nTrials = ir.mosaic{1}.numberTrials;
        for itrial = 1:nTrials
            fprintf('Trial %d\n',itrial);
            ir = irComputeSpikes(ir);
        end
end

end



