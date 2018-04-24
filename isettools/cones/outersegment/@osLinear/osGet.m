function val = osGet(obj, varargin)
% Gets isetbio outersegment object parameters.
% 
% Syntax:
%   val = osGet(obj, [varargin])
%
% Description:
%    Get the isetbio outersegment object parameter by name.
%
%    Examples are contained in the code. To access, type 'edit osGet.m'
%    into the Command Window.
%
% Inputs:
%	 obj      - The outersegment object to retrieve parameter(s) from.
%    varargin - The parameter and additional information as required. The
%               list of possible parameters to retrieve includes:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'sConeFilter'} - the linear filter for S-cone temporal response
%       {'mConeFilter'} - the linear filter for M-cone temporal response
%       {'lConeFilter'} - the linear filter for L-cone temporal response
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
%
% Outputs:
%    val      - The value of the requested parameter.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Cone current signal does not appear to be supported any
%      longer? Should we remove this option? And the size option since it
%      references conecurrentsignal? Also the conefilters (l, m, and s)?]
%

% History:
%    08/xx/15  JRG NC DHB  Created
%    02/13/18  jnm         Formatting

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

%% Make key properties that can be set required arguments.
allowableFieldsToSet = {...
    'noiseflag', ...
    'sconefilter', ...
    'mconefilter', ...
    'lconefilter', ...
    'patchsize', ...
    'timestep', ...
    'size', ...
    'conecurrentsignal'};
p.addRequired('what', @(x) any(validatestring(ieParamFormat(x), ...
    allowableFieldsToSet)));

%% Parse and put results into structure p.
p.parse(varargin{:}); 
params = p.Results;

switch ieParamFormat(params.what)
    case{'sconefilter'}
        val = obj.sConeFilter;
        
    case{'mconefilter'}
        val = obj.mConeFilter;
    
    case{'lconefilter'}
        val = obj.lConeFilter;        
        
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'size'}
        val = size(obj.coneCurrentSignal);
        
    case{'noiseflag'}
        val = obj.noiseFlag;
        
    case{'conecurrentsignal'}
        val = obj.coneCurrentSignal;
end
