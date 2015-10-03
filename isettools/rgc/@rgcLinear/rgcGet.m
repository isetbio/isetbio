function val = rgcGet(obj, varargin)
% rgcGet: a method of @rgcLinear that gets rgcLinear object 
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
        'animal',...
        'input',...
        'eyeLeftOrRight',...
        'patchLocationPolarRadiusMicrometers',...
        'patchLocationPolarAngleDegrees',...
        'temporalEquivEcc',...        
        'numberCellTypes',...
        'namesCellTypes',...        
        'mosaic',...        
        'noiseFlag'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));

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
    
    case{'input'}
        val = obj.input;
    case{'animal'}
        val = obj.animal;
    case{'eyeLeftOrRight'}        
        val = obj.eyeLeftOrRight;
    case{'patchLocationPolarRadiusMicrometers'}        
        val = obj.patchLocationPolarRadiusMicrometers;
    case{'patchLocationPolarAngleDegrees'}        
        val = obj.patchLocationPolarAngleDegrees;
    case{'temporalEquivEcc'}        
        val = obj.temporalEquivEcc;
    case{'numberCellTypes'}
        val = obj.numberCellTypes;
    case{'namesCellTypes'}        
        val = obj.namesCellTypes;
    case{'mosaic'}        
        val = obj.mosaic;        
    case{'noiseFlag'}            
        val = obj.noiseFlag;
        
end

