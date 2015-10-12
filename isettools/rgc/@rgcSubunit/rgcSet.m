function obj = rgcSet(obj, varargin)
% rgcSet: a method of @rgcLNP that sets rgcLNP object 
% parameters using the input parser structure.
% 
% Parameters:
%       {''} -  
% 
% 9/2015 JRG 

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 

% % % We could do set using the superclass method
% obj = mosaicSet@rgcMosaic(obj, varargin{:});

% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...  
        'name',...
        'input',...
        'temporalEquivEcc',...       
        'mosaic'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
p.addRequired('value');

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch lower(params.what)   
 
    case{'name'}
        obj.name = params.value;
    case{'input'}
        obj.input = params.value;
    case{'temporalequivecc'}        
        obj.temporalEquivEcc = params.value;
    case{'mosaic'}        
        obj.mosaic = params.value;      
end

