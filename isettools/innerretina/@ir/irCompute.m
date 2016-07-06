function ir = irCompute(ir, input, varargin)
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
% JRG (c) isetbio team

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
vFunc = @(x) ~isempty(validatestring(class(x),{'osDisplayRGB','bipolar'}));
p.addRequired('input',vFunc);

p.parse(ir,input,varargin{:});

%% Linear stage of the computation
ir = irComputeLinearSTSeparable(ir, input);
% irPlot(ir,'response linear');

%% Compute spikes for each trial
switch class(ir.mosaic{1})
    case {'rgcLinear'};
        % No linear response implemented yet.
        disp('No spikes computed for linear RGC mosaic');   
    otherwise
        % Runs for rgcLNP, rgcGLM
        ir = irComputeSpikes(ir);
end

end