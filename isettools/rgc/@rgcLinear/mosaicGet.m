function val = mosaicGet(obj, param, varargin)
% @rgcLinear subclass mosaicGet, superclass is @rgcMosaic 
%
% parameters using the input parser structure.
% 
%       val = mosaicGet(rgc.mosaic, param, property)
% 
% Inputs: rgc object, property to be gotten
% 
% Outputs: val of property
% 
% Properties special to the rgcLinear class
%   N/A yet
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'linearResponse')
% 
% 9/2015 JRG 

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% The permissible fields for this class are below.  In addition, there are
% parameters in the superclass, rgcMosaic.
p.addRequired('param',@ischar);

% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = p.Results.param;

%% Set key-value pairs.
switch ieParamFormat(param)
    
    % The special cases for rgcLinear should be here
    
    otherwise
        val = mosaicGet@rgcMosaic(obj,param,varargin{:});
                 
end

end
