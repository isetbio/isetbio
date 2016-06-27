function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters.
% 
% Parameters:
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'photonRate'} - photon rate from sensor, copied for osIdentity.
% 
% noiseFlag = 0;
% adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);
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
allowableFieldsToSet = {'noiseflag','photonrate','patchsize','timestep'};
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addRequired('value');

% Define what units are allowable.
allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps

% Set up key value pairs.
% Defaults units:
p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what);  % Lower case and remove spaces

   
    case{'patchsize'}
        obj.patchSize = params.value;
        
    case{'timestep'}
        obj.timeStep = params.value;
        
    case{'photonrate'}
        obj.rgbData = params.value;
               

end

