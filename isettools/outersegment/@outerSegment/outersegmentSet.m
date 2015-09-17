function obj = outersegmentSet(obj, varargin)
% outersegmentSet: a method of @outersegment that sets  outersegment object 
% parameters using the input parser structure.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
% 
% 
% 8/2015 JRG NC DHB


% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {'noiseflag'};
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
p.addRequired('value');

% Define what units are allowable.
allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps

% Set up key value pairs.
% Defaults units:
p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch lower(params.what);  % Lower case and remove spaces

    
    case{'noiseflag'}
        obj.noiseFlag = params.value;
               

end

