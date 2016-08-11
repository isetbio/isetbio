function val = mosaicGet(obj, param, varargin)
% Gets a property from an rgcLinear object.
% 
%       val = mosaicGet(rgc.mosaic, param, property)
% 
% The rgcLinear mosaic object has all of its properties defined by the
% rgcMosaic superclass. A call to mosaicGet for the rgcLinear object will
% pass the mosaicGet call onto the superclass method.
% 
% Inputs: 
%   obj    - rgc object
%   param  - parameter string
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%   val - parameter value
% 
%  Properties that can be gotten: see @rgcMosaic/mosaicGet.m
% 
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'linearResponse')
% 
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% The permissible fields for this class are below.  In addition, there are
% parameters in the superclass, rgcMosaic.
p.addRequired('param',@ischar);

p.parse(param,varargin{:}); 
param = p.Results.param;

%% Set key-value pairs.
switch ieParamFormat(param)
    
    % All properties for rgcLinear are defined in the rgcMosaic superclass,
    % so see @rgcMosaic/mosaicGet.m.    
    otherwise
        val = mosaicGet@rgcMosaic(obj,param,varargin{:});
                 
end

end
