function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters.
% 
% Syntax:
%   obj = osSet(obj, varargin)
%
% Description:
%    Set the isetbio outer segment object properties for the base class.
%
%    Examples are contained in the code. To access, type 'edit osSet.m'
%    into the Command Window.
%
% Inputs:
%    obj      - The outer segment object you wish to assign/change a
%               property value of.
%
% Outputs:
%    obj      - The modified outer segment object.
%
% Optional key/value pairs:
%    varargin - contains a key/value pair of the parameter and it's
%               associated value that you wish to set. The possible options
%               for the parameter include:
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'photonRate'} - photon rate from sensor, copied for osIdentity.
%       {'noiseFlag'} - 'random'-Default, 'frozen', 'none'.
%
% Notes:
%    * [Note: JNM - Changed the listed value of noiseFlag to an appropriate
%      value to allow the example to work. Added noiseFlag case below, for
%      value to actually set!]
%

% History:
%    08/xx/15  JRG NC DHB  Created
%    02/14/18  jnm         Formatting

% Examples:
%{
    adaptedOS = osCreate('identity');
    noiseFlag = 'frozen';
    adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag)
%}


% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag', ...
    'photonrate', ...
    'patchsize', ...
    'timestep'};
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:});
params = p.Results;

switch ieParamFormat(params.what)  % Lower case and remove spaces
    case{'noiseflag'}
        obj.noiseFlag = params.value;
        
    case{'patchsize'}
        obj.patchSize = params.value;
        
    case{'timestep'}
        obj.timeStep = params.value;
        
    case{'photonrate'}
        obj.photonRate = params.value;
               
end
