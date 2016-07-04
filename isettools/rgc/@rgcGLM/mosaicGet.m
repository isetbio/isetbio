function val = mosaicGet(obj, param, varargin)
% mosaicGet for @rgcGLM, a subclass of @rgcMosaic
%
%   val = mosaicGet(rgc.mosaic, param, property)
% 
% Inputs: rgc object, property to be gotten
% 
% Outputs: 
%    val of param
% 
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'psthResponse')
% 
% ISETBIO Team, 2016

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
    
    case{'generatorfunction'}
        val = obj.generatorFunction;
    case{'numbertrials'}
        % Is this in the super class?
        %  val = obj.numberTrials;
        val = size(obj.responseSpikes,3);
    case{'responsevoltage'}
        val = obj.responseVoltage;
    case{'postspikefilter'}
        val = obj.postSpikeFilter;
    case{'couplingfilter'}
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        val = obj.couplingMatrix;
        
    otherwise
        % Super class
        val = mosaicGet@rgcMosaic(obj,param,varargin{:});
end

end

