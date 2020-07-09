function val = eyeGet(obj,param,varargin)
% The Get function for sceneEye parameters
%
% Synopsis
%   val = eyeGet(obj,param,varargin)
% 
% Brief description:
%   These gets used to be in sceneEye, but as we integrated with iset3d we
%   moved them out to here
% 
% Input
%  obj:     The sceneEye object with obj.recipe slot
%  param:   Name of the parameter (case and spaces ignored)
%  
% Optional key/val pairs
%  N/A
%
% Output
%  val:  Returned value
%
% See also
%   eyeSet, recipeGet, recipeSet
%

switch ieParamFormat(param)
    case 'angularsupport'
        ss = obj.sampleSize;
        chordSpatialSamples = (1:obj.resolution).*ss - obj.width/2;
        val = atand(chordSpatialSamples/obj.distance2chord);
    case 'samplesize'
        val = obj.width / obj.resolution;
        
    case 'width'
        val = obj.width;
        
    case 'height'
        val = obj.width;
        
    case 'distance2chord'
        val = 2 * tand(obj.fov / 2) * (obj.distance2chord);
        
    otherwise
        % Pass through to the main get of the recipe.
        val = obj.recipe.get(param,varargin{:});
end

end
        
