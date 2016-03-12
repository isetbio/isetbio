function val = osGet(obj, param)
% Get base class outer segment parameters
% 
% Parameters:
%  {'noiseFlag'} -  gets noise flag, noise-free ('0') or noisy ('1')
%  {'ConeCurrentSignal'} - cone current as a function of time
% 
% 7/2015 JRG NC DHB

if ~exist('param','var'), error('Parameter required'); end

switch ieParamFormat(param);  % Lower case and remove spaces

    case {'noiseflag'}        
        val = obj.noiseFlag;
        
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'conecurrentsignal'}
        val = obj.coneCurrentSignal;
        
    otherwise
        error('Unknown parameter %s\n',param);
       
end

