function obj = osSet(obj, varargin)
% Sets isetbio outersegment object parameters
%
% Multiple parameter/value pairs can be set at a time
%
% Parameters:
%   sConeFilter - the linear filter for S-cone temporal response
%   mConeFilter - the linear filter for M-cone temporal response
%   lConeFilter - the linear filter for L-cone temporal response
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
        case {'sconefilter'}
            % Temporal impulse response
            obj.sConeFilter = value;
            
        case {'mconefilter'}
            % Temporal impulse response
            obj.mConeFilter = value;
            
        case {'lconefilter'}
            % Temporal impulse response
            obj.lConeFilter = value;            
            
        otherwise
            % If not part of this class, check, the parent class.
            obj = osSet@outerSegment(obj,param,value);
    end
    
end

end