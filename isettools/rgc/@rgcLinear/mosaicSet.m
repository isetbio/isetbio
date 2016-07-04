function obj = mosaicSet(obj, param, val, varargin)
% rgcLinear subclass mosaic set.  Superclass is @rgcMosaic
% 
%   rgc.mosaic = mosaicSet(rgc.mosaic, property, value, varargin)
%  
% Inputs: rgc object, 
%   property value pair
% 
% Outputs: 
%  obj - with property set appropriately
% 
% Examples:
%    mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%    mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% 9/2015 JRG 

%% Parse
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

p.addRequired('param',@ischar);
p.addRequired('val');

p.parse(param, val, varargin{:}); 
param = p.Results.param;
val   = p.Results.val;

%% Set key-value pairs.

% @JRG: Please add comments about what these parameters are and their
% potential values.
switch ieParamFormat(param)

    % Special to this class should be here
    
    otherwise
        mosaicSet@rgcMosaic(obj,param,val,varargin{:});

        % DELETE ME
%     case{'celltype'}
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
 
end

