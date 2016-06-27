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

%% Loop through param/value pairs

for ii=1:2:length(varargin)

    param = ieParamFormat(varargin{ii});
    value = varargin{ii+1};
    
    switch param
                    
        case {'noiseflag'}
            obj.noiseFlag = value;         
                      
        case{'patchsize'}
            obj.patchSize = params.value;
            
        case{'timestep'}
            obj.timeStep = params.value;            
        
        case{'conecurrentsignal'}
            obj.coneCurrentSignal = value;
            
    end
    
end

