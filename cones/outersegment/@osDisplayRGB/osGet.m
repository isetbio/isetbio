function val = osGet(obj, varargin)
% Gets the isetbio outersegment object parameters.
% 
% Syntax:
%    val = osGet(obj, varargin)
%
% Description:
%    Retrieve the base class outer segment parameters
%
%    Examples are contained in the code. To access, type 'edit osGet.m'
%    into the Command Window.
%
% Inputs:
%    obj      - The outer segment object
%    varargin - The parameter you wish to retrieve the value of. Options
%               include the following:
%        'rgbData'    - Scene RGB data to pass to "black box" RGC GLM Model
%        'timestep'   - Noisy cone current signal
%        'patchSize'  - Cone current as a function of time
%        'size'       - Array size of RGB Data
%
% Outputs:
%    val   - The value of the requested parameter
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/xx/15  JRG  Created
%    02/14/18  jnm  Formatting

% Examples:
%{
    adaptedOS = osCreate;
   osGet(adaptedOS, 'noiseFlag')
%}

%% Check for the number of arguments and create parser object.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% Make key properties that can be set required arguments.
allowableFieldsToSet = {...
    'noiseflag', ...
    'rgbdata', ...
    'patchsize', ...
    'timestep', ...
    'size'};
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:});
params = p.Results;

switch ieParamFormat(params.what)  % Lower case and remove spaces
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'rgbdata'}
        val = obj.rgbData;
        
    case{'size'}
        val = size(obj.rgbData);

end
