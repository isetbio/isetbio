function val = osGet(obj, varargin)
% Gets the isetbio outersegment object parameters.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% osGet(adaptedOS, 'noiseFlag')
% 
% 8/2015 JRG NC DHB


%% Check for the number of arguments and create parser object.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments.
allowableFieldsToSet = {...
    'noiseflag',...
    'patchsize',...
    'timestep',...
    'size',...
    'conecurrentsignal',...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(params.what);  % Lower case and remove spaces

    case {'noiseflag'}        
        val = obj.noiseFlag;
        
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'size'}
        val = size(obj.coneCurrentSignal);
        
    case{'conecurrentsignal'}
        val = obj.coneCurrentSignal;
end

