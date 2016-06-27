function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters.
%
% Parameters:
% 
%       {'sConeFilter'} - the linear filter for S-cone temporal response
%       {'mConeFilter'} - the linear filter for M-cone temporal response
%       {'lConeFilter'} - the linear filter for L-cone temporal response
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% Example:
%
%   adaptedOS = osSet(adaptedOS, 'noiseFlag', 0, 'time step',0.001);
%
% 8/2015 JRG NC DHB


%% Set all the key-value pairs.

for ii=1:2:length(varargin)
    param = ieParamFormat(varargin{ii});
    value = varargin{ii+1};
    switch param        
        
        case {'noiseflag'}
            obj.noiseFlag = value;
        
        case {'sconefilter'}
            % Temporal impulse response
            obj.sConeFilter = value;
            
        case {'mconefilter'}
            % Temporal impulse response
            obj.mConeFilter = value;
            
        case {'lconefilter'}
            % Temporal impulse response
            obj.lConeFilter = value;            
                      
        case{'patchsize'}
            obj.patchSize = params.value;
            
        case{'timestep'}
            obj.timeStep = params.value;            
        
        case{'conecurrentsignal'}
            obj.coneCurrentSignal = value;
    end
    
end

end