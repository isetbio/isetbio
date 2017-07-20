function obj = set(obj, param, val, varargin)
%  Sets a property for an rgcGLM object.
%
%   @rgcGLM.set(param, val, varargin)
%
% Inputs: 
% 
%   obj    - rgc object
%   param  - parameter string
%   val    - parameter value
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%    obj with property set appropriately
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated
%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

p.addRequired('param');
p.addRequired('val');

p.parse(param,val,varargin{:});
param = ieParamFormat(p.Results.param);
val   = p.Results.val;

switch param
    
    % GLM subclass parameters
    case{'generatorfunction'}
        % An inline function, e.g. @exp
        obj.generatorFunction = val;
    case{'numbertrials'}
        % The number of trials for which spikes are computed from the
        % linear response
        obj.numberTrials = val;
    case{'responsevoltage'}        
        % The "membrane voltage" from the spike computation in units of
        % conditional intensity. This signal contains the effects of the
        % post spike filter and coupling filters for individual spikes in a
        % given trial.
        obj.responseVoltage = val;
    case{'postspikefilter'}
        % The post spike filter, in units of conditional intensity
        obj.postSpikeFilter = val;
    case{'couplingfilter'}
        % The coupling filter in units of conditional intensity.
        obj.couplingFilter = val;
    case{'couplingmatrix'}
        % The weights on the coupling filter between cells in normailzed
        % units between zero and one.
        obj.couplingMatrix = val;
    otherwise
        % Superclass
        set@rgcMosaic(obj,param,val,varargin{:});
        
end
end