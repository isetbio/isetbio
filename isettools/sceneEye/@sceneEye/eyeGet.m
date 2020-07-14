function val = eyeGet(obj,param,varargin)
% The Get function for sceneEye and ISET3d recipe parameters
%
% Synopsis
%   val = eyeGet(obj,param,varargin)
% 
% Brief description:
%   Get parameters from the sceneEye class or from the ISET3ed recipe
%   attached to the sceneEye object.
% 
% Input
%  obj:     The sceneEye object with obj.recipe defining the rendering
%  param:   Name of the parameter (case and spaces ignored).  The param can
%           be from the sceneEye class or the ISET3d recipe class.
%  
% Optional key/val pairs
%  N/A
%
% Output
%  val:  Returned value
%
% Parameters
%
%    'name'  - Name of this sceneEye object
%    'model name' - Name of the lens model (default: navarro)
%    'recipe'     - ISET3d recipe
%    'use pinhole' ('use optics')  - Render a scene (use pinhole) or
%                  optical image ('use optics')
%    'fov'  - Sets the  horizontal field of view (deg) by adjusting the
%             retinaSemidiam of the retinal geometry.
%    'lens density' - Lens pigment density
%
% See also
%   eyeSet, recipeGet, recipeSet
%

%%
switch ieParamFormat(param)
    case 'name'
        % thisEye.get('name')
        % Name of this object
        val = obj.name;
        
    case 'modelname'
        % thisEye.get('model name')
        % Name of the eye model
        val = obj.modelname;
        
    case 'angularsupport'
        warning('angular support - Not likely to work')
        ss = obj.sampleSize;
        chordSpatialSamples = (1:obj.resolution).*ss - obj.width/2;
        val = atand(chordSpatialSamples/obj.distance2chord);
        
    case 'usepinhole'
        % This parameter forces the rendering to swap out the camera and
        % use pinhole optics instead.  It is used for debugging and testing
        % because the rendering can be pretty fast.
        val = obj.usePinhole;
        
    case 'useoptics'
        % This parameter forces the rendering to swap out the camera and
        % use pinhole optics instead.  It is used for debugging and testing
        % because the rendering can be pretty fast.
        val = ~obj.usePinhole;
        
    case 'lensdensity'
        % Adjust the lens pigment density
        val = obj.lensDensity;
        
    case 'recipe'
        % The rendering recipe
        val = obj.recipe;
        
    otherwise
        % Pass through to the get of the ISET3d recipe.
        val = obj.recipe.get(param,varargin{:});
end

end
        
