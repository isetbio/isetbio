function obj = osSet(obj, varargin)
% Sets isetbio outersegment object for biophys model
%
%
%
% Example:
%   adaptedOS = osSet(adaptedOS, 'noise flag', 0);
% 
% 8/2015 JRG NC DHB

%% Loop through param/value pairs

for ii=1:2:length(varargin)

    param = ieParamFormat(varargin{ii});
    value = varargin{ii+1};
    
    switch param
        
        otherwise
            % If not part of this class, check, the parent class.
            obj = osSet@outerSegment(obj,param,value);
    end
    
end

