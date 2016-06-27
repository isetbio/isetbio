function obj = osSet(obj, varargin)
% Sets the isetbio outersegment object parameters.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% osGet(adaptedOS, 'noiseFlag')

%%
% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'noiseflag',...
    'conecurrentsignal',...
    'patchsize',...
    'timestep'};
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addRequired('value');

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;
%% Loop through param/value pairs

switch ieParamFormat(params.what);  % Lower case and remove space
                    
        case {'noiseflag'}
            obj.noiseFlag = params.value;         
                      
        case{'patchsize'}
            obj.patchSize = params.value;
            
        case{'timestep'}
            obj.timeStep = params.value;            
        
        case{'conecurrentsignal'}
            obj.coneCurrentSignal = params.value;
 
end

