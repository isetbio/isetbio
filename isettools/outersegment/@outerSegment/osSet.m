function obj = osSet(obj, param, value)
% Sets isetbio outersegment object properties for base class
% 
% Parameters:
%   {'noise flag'} -  sets current as noise-free ('0') or noisy ('1')
%   {'time step'}  -  
%   {'patch size'} -  
%   {'cone current signal'} -  
%
% Example:
%   adaptedOS = osSet(adaptedOS, 'noise flag', 0);
% 
% 8/2015 JRG NC DHB

%% Check for the number of arguments and create parser object.
if ~exist('param','var'), error('Parameter required'); end
if ~exist('value','var'), error('Value required'); end

%%
switch ieParamFormat(param)
    
    case{'noiseflag'}
        obj.noiseFlag = value;               
               
    case{'timestep'}
        obj.timeStep = value;
        
    case{'patchsize'}
        % Spatial sample spacing
        obj.patchSize = value;
        
    case{'conecurrentsignal'}
        obj.coneCurrentSignal = value; 
end
