function val = eyeGet(obj,param,varargin)
% The Get function for sceneEye parameters
%
% Synopsis
%   val = eyeGet(obj,param,varargin)
% 
% Brief description:
%   Get parameters from the sceneEye class or from the recipe attached to
%   the sceneEye object.
% 
% Input
%  obj:     The sceneEye object with obj.recipe defining the rendering
%  param:   Name of the parameter (case and spaces ignored).  The param can
%           be from the sceneEye class or the recipe class.
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

%%
switch ieParamFormat(param)
    case 'angularsupport'
        ss = obj.sampleSize;
        chordSpatialSamples = (1:obj.resolution).*ss - obj.width/2;
        val = atand(chordSpatialSamples/obj.distance2chord);
        
    case 'width'
        val = obj.width;
        
    case 'height'
        val = obj.width;
   
    case 'usepinhole'
        % This parameter forces the rendering to swap out the camera and
        % use pinhole optics instead.  It is used for debugging and testing
        % because the rendering can be pretty fast.
        val = obj.usePinhole;
        
    case 'recipe'
        % The rendering recipe
        val = obj.recipe;
        
    otherwise
        % Pass through to the main get of the recipe.
        val = obj.recipe.get(param,varargin{:});
end

end
        
