function [ir, nTrialsSpikes] = irCompute(ir, input, varargin)
% Computes the rgc mosaic responses to an input
%
%    ir = irCompute(ir, input, varargin)
%
% Inputs:
%   ir: inner retina object
%   input - There are two types of possible input objects.
%
%      'osDisplayRGB' - frame buffer values for a spatiotemporal stimulus
%           stored in an outer segment object.
%      'bipolar' - the bipolar cell object with a signal that has gone
%           through temporal filtering and possibly spatial subunit
%           nonlinearities.
%
% Computes the linear and spike responses for all the mosaics within the
% inner retina object.  (Except: There are not spikes for the rgcLinear
% class).
% 
% For each mosaic, a space-time separable linear response is computed. This
% stage of the computation is stored in 'responseLinear'.  This is managed
% in irComputeLinearSTSeparable.  There is no noise added in the linear
% part.
%
% The spikes are computed irComputeSpikes routine. The spiking can have a
% random element.  So, we may run the conversion from linear to spikes
% multiple times, effectively producing spike rasters.
%
% Outputs:
%  ir: the inner retina object with responses attached to each mosaic
%
% Example:
%   ir.compute(identityOS);
%   irCompute(ir,identityOS);
%
% See also: rgcMosaic, irComputeSpikes, irComputeLinearSTSeparable
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated
%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
if length(input) == 1
    vFunc = @(x) ~isempty(validatestring(class(x),{'osDisplayRGB','bipolar'}));
else
    vFunc = @(x) ~isempty(validatestring(class(x{1}),{'osDisplayRGB','bipolar'}));
end
p.addRequired('input',vFunc);
p.addParameter('coupling',false,@islogical);
p.addParameter('nTrialsInput',  [], @isnumeric);

p.parse(ir,input,varargin{:});
coupling = p.Results.coupling;

nTrialsInput = p.Results.nTrialsInput; 

%% Linear stage of the computation

if ~isempty(nTrialsInput)  
    [ir,nTrialsLinearResponse] = irComputeLinearSTSeparable(ir, input, 'nTrialsInput',nTrialsInput);
else
    ir = irComputeLinearSTSeparable(ir, input);
end
% irPlot(ir,'response linear');

%% Compute spikes for each trial
switch class(ir.mosaic{1})
    case {'rgcLinear'};
        % No linear response implemented yet.
        disp('No spikes computed for linear RGC mosaic');   
    otherwise
        % Runs for rgcLNP, rgcGLM
        % Check the coupling field to decide on the coupling parameter
        if ~isempty(nTrialsInput) 
            [ir, nTrialsSpikes] = irComputeSpikes(ir,'coupling',coupling,'nTrialsLinearResponse',nTrialsLinearResponse);
        else
            [ir, nTrialsSpikes] = irComputeSpikes(ir,'coupling',coupling);
        end
end

end