function [ir, nTrialsSpikes] = irCompute(ir, bp, varargin)
% IRCOMPUTE - Computes the rgc mosaic responses to an input
%
%   ir = irCompute(ir, bipolar, varargin)
%
% Required inputs
%  'ir' - inner retina object
%  'bipolar' - the bipolar mosaic object
%
% Optional inputs
%   'bipolarTrials' -
%
% Computes the continuous (linear) and spike responses for each of the
% mosaics within the inner retina object.  (Note: There are no spikes
% for the rgcLinear class).
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
%  nTrialsSpikes:  (trials x xPos x yPos x Time)
%     Binary matrix indicating the spike times for all the trials.
%     The last one is stored in the mosaics of the inner retina object in
%     the responseSpikes slot.
%
% Example:
%   ir.compute(bipolar);
%   irCompute(ir,bipolar);
%
% See also: rgcMosaic, irComputeSpikes, irComputeLinearSTSeparable
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated
%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
vFunc = @(x) (isequal(class(x),'bipolar')||isequal(class(x{1}),'bipolar'));
p.addRequired('bp',vFunc);

p.addParameter('coupling',false,@islogical);
p.addParameter('bipolarTrials',  [], @isnumeric);  % Multiple bipolar trials

p.parse(ir,bp,varargin{:});
coupling = p.Results.coupling;

bipolarTrials = p.Results.bipolarTrials; 

%% Linear stage of the computation

if ~isempty(bipolarTrials)  
    [ir,nTrialsLinearResponse] = irComputeLinearSTSeparable(ir, bp, 'bipolarTrials',bipolarTrials);
else
    ir = irComputeLinearSTSeparable(ir, bp);
end
% irPlot(ir,'response linear');

%% Compute spikes from linear response; do for multiple trials

% This should be for ii=1:length(ir.mosaic)
switch class(ir.mosaic{1})
    case {'rgcLinear'}
        % No linear response implemented yet.
        disp('No spikes computed for linear RGC mosaic');   
    otherwise
        % Runs for rgcLNP, rgcGLM
        % Check the coupling field to decide on the coupling parameter
        if ~isempty(bipolarTrials) 
            [ir, nTrialsSpikes] = irComputeSpikes(ir,'coupling',coupling,'nTrialsLinearResponse',nTrialsLinearResponse);
        else
            [ir, nTrialsSpikes] = irComputeSpikes(ir,'coupling',coupling);
        end
end

end