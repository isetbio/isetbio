function val = mosaicGet(obj, param, varargin)
% mosaicGet for LNP subclass; superclass is @rgcMosaic
%
%   val = mosaicGet(rgc.mosaic, param, varargin)
%
% Inputs: 
%   rgc object
%   param - to retrieve
%   vararing depends on parameter
%
% Outputs: 
%  val of parameter
%
% LNP parameters
%
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cell type')
%   val = mosaicGet(rgc1.mosaic{3}, 'psth response')
%

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;
p.KeepUnmatched = true;
p.addRequired('param');

% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = p.Results.param;

%% Set key-value pairs.
switch ieParamFormat(param)
    
    % Specific to the GLM case
    case{'generatorfunction'}
        val = obj.generatorFunction;
    case{'postspikefilter'}
        val = obj.postSpikeFilter;
    case{'numbertrials'}
        if ~isempty(obj.responseSpikes)
            val = size(obj.responseSpikes,3);
        else
            val = 0;
        end
    case{'responsevoltage'}
        val = obj.responseVoltage;
    case{'couplingfilter'}
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        val = obj.couplingMatrix;
        
    otherwise
        val = mosaicGet@rgcMosaic(obj,param,varargin{:});
end      

end

