function ir = irCompute(ir, input, varargin)
% Computes the rgc mosaic responses to an input
%
%
% The input can be of several types, for historical reasons.
%
% Inputs:
%   ir: inner retina object,
%   input - one of
%      {'osDisplayRGB','osIdentity','osLinear','osBioPhys','bipolar'}
%
% Computes the responses for all the mosaics attached to the inner retina
% object.
%
% For each mosaic, the linear response is computed by spatial convolution
% of the center and surround RFs. Then, the temporal responses for the
% center and surround is computed, separably. This stage of the computation
% is stored in responseLinear.
%
% The linear response is sent to the irComputeSpikes routine.  Each rgc
% type will have its own compute method.  For all types the initial linear
% method is carried out in here.
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
% We plan to move the linear computation to a separate routine,
% irComputeContinuous (or Linear) and carry out the first part of this code
% in that routine, rather than as the switch statement that is here now
%

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
p.addRequired('input',@(x) ~isempty(validatestring(class(x),{'osDisplayRGB','osIdentity','osLinear','osBioPhys','bipolar'})));

p.parse(ir,input,varargin{:});

% Linear stage of the computation
ir = irComputeContinuous(ir, input);

% Compute spikes for each trial
switch class(ir.mosaic{1})
    case {'rgcLinear'};
        % No nonlinear response
        error('Not yet implemented');
    case{'rgcPhys'}
        
    otherwise
        nTrials = ir.mosaic{1}.numberTrials;
        for itrial = 1:nTrials
            itrial
            ir = irComputeSpikes(ir);
        end
end



