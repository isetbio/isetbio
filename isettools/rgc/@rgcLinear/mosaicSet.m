function obj = mosaicSet(obj, param, val, varargin)
%  Sets a property for an rgcLinear object.
% 
%   rgc.mosaic = mosaicSet(rgc.mosaic, property, value, varargin)
% 
% The rgcLinear mosaic object has all of its properties defined by the
% rgcMosaic superclass. A call to mosaicSet for the rgcLinear object will
% pass the mosaicSet call onto the superclass method.
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
%    mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%    mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;
p.addRequired('param',@ischar);
p.addRequired('val');

p.parse(param, val, varargin{:}); 
param = p.Results.param;
val   = p.Results.val;

%% Set key-value pairs.

switch ieParamFormat(param)

    % All properties for rgcLinear are defined in the rgcMosaic superclass,
    % so see @rgcMosaic/mosaicSet.m.    
    otherwise
        mosaicSet@rgcMosaic(obj,param,val,varargin{:});
 
end

