function obj = mosaicSet(obj, param, val, varargin)
% rgcLinear subclass mosaic.  Superclass is @rgcMosaic
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
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

p.addRequired('param',@ischar);
p.addRequired('val');

p.parse(param, val, varargin{:}); 
param = p.Results.param;
val   = p.Results.val;

%% Set key-value pairs.

switch ieParamFormat(param)

    % Special to this class should be here
    
    otherwise
        mosaicSet@rgcMosaic(obj,param,val,varargin{:});
 
end

