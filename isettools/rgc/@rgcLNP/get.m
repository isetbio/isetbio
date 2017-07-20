function val = get(obj, param, varargin)
% Gets a rgcLNP parameter
%
%   val = @rgcLNP.get(param, varargin)
%
% Passes along parameters not found here to @rgcMosaic.get
%
% Inputs: 
%   obj      - rgc object
%   param    - parameter string
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
%   val = @rgcLNP.get('cell type')
%   val = @rgcLNP.get('psth response')
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;
p.addRequired('param');

% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = p.Results.param;

%% Set key-value pairs.
switch ieParamFormat(param)
    
    % Specific to the LNP case
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
        
    otherwise
        val = get@rgcMosaic(obj,param,varargin{:});
end      

end

