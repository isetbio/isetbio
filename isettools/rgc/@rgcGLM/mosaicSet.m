function obj = mosaicSet(obj, param, val, varargin)
% rgcMosaicSet: a method of @rgcMosaic that sets rgcMosaic object
% parameters using the input parser structure.
%
%       rgc.mosaic = mosaicSet(rgc.mosaic, param, val, varargin)
%
% Inputs: rgc object, key-value pair of property and value to which it is
%   being set.
%
% Outputs: obj with property set appropriately.
%
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
%
% 9/2015 JRG

%% Parse
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

p.addRequired('param');
p.addRequired('val');

p.parse(param,val,varargin{:});
param = ieParamFormat(p.Results.param);
val   = p.Results.val;

% Set key-value pairs.
switch param
    
    % GLM specific
    case{'generatorfunction'}
        obj.generatorFunction = val;
    case{'numbertrials'}
        obj.numberTrials = val;
    case{'responsevoltage'}
        obj.responseVoltage = val;
    case{'postspikefilter'}
        obj.postSpikeFilter = val;
    case{'couplingfilter'}
        obj.couplingFilter = val;
    case{'couplingmatrix'}
        obj.couplingMatrix = val;
    otherwise
        % Superclass
        mosaicSet@rgcMosaic(obj,param,val,varargin{:});
        
end
end


% Delete me
%  case{'celltype'}
%         obj.cellType = val;
%     case{'rfdiameter'}
%         obj.rfDiameter = val;
%     case{'rfdiamagnitude'}
%         obj.rfDiaMagnitude = val;
%     case{'celllocation'}
%         obj.cellLocation = val;
%     case{'srfcenter'}
%         obj.sRFcenter = val;
%     case{'srfsurround'}
%         obj.sRFsurround = val;
%     case{'tcenter'}
%         obj.tCenter = val;
%     case{'tsurround'}
%         obj.tSurround = val;
%     case{'responselinear'}
%         obj.responseLinear = val;
%
