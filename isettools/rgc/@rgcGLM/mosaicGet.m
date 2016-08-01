function val = mosaicGet(obj, param, varargin)
% Gets a property from an rgcGLM object.
%
%   val = mosaicGet(rgc.mosaic, param, varargin)
%
% Inputs: 
%   obj    - rgc object
%   param  - parameter string
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%   val - parameter value
% 
%  Properties that can be gotten: 
%     generatorFunction
%     postSpikeFilter
%     numberTrials
%     responseVoltage
%     couplingFilter
%     couplingMatrix
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cell type')
%   val = mosaicGet(rgc1.mosaic{3}, 'psth response')
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

p.addRequired('param');
p.parse(param,varargin{:});
param = ieParamFormat(p.Results.param);

%% Set key-value pairs.
switch param
    % Specific to the GLM case
    case{'generatorfunction'}
        % An inline function, e.g. @exp
        val = obj.generatorFunction;
    case{'postspikefilter'}
        % The post spike filter, in units of conditional intensity
        val = obj.postSpikeFilter;
    case{'numbertrials'}
        % The number of trials for which spikes are computed from the
        % linear response
        if ~isempty(obj.responseSpikes)
            val = size(obj.responseSpikes,3);
        else
            val = 0;
        end
    case{'responsevoltage'}
        % The "membrane voltage" from the spike computation in units of
        % conditional intensity. This signal contains the effects of the
        % post spike filter and coupling filters for individual spikes in a
        % given trial.
        val = obj.responseVoltage;
    case{'couplingfilter'}
        % The coupling filter in units of conditional intensity.
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        % The weights on the coupling filter between cells in normailzed
        % units between zero and one.
        val = obj.couplingMatrix;
        
    otherwise
        % Super class
        val = mosaicGet@rgcMosaic(obj,param,varargin{:});
end

end

