function obj = eyeSet(obj,param,val,varargin)
% Set methods for sceneEye and ISET3d recipe dependent variables
%
% Synopsis
%   obj = eyeSet(obj,param,val,varargin)
%
% Brief description
%   Set the sceneEye or iset3d recipe parameters through this call.  Case
%   and spaces are ignored in the param string.
%
% Input
%  obj:   sceneEye class
%  param: string defining the parameter
%  val:   Value of the parameter
%
% 
% Optional key/val parmeters
%  N/A
%
% Output
%   obj:  Modified sceneEye class
%
% Settable parameters in eyeSet
%    'name'  - Name of this sceneEye object
%    'model name' - Name of the lens model (default: navarro)
%    'recipe'     - ISET3d recipe
%    'use pinhole' ('use optics')  - Render a scene (use pinhole) or
%                  optical image ('use optics')
%    'fov'  - Sets the  horizontal field of view (deg) by adjusting the
%             retinaSemidiam of the retinal geometry.
%    'lens density' - Lens pigment density
%s
%   Many more parameters can be set in the recipe.  See recipeSet.
%
% See also
%  sceneEye, recipeSet/Get

%%  Force param to lower case, no spaces\
param = ieParamFormat(param);

%%  Main switch statment
switch param
    case 'name'
        obj.name = val;

    case 'modelname'
        % If we set the eye model, we need to change the retina distance and
        % radius.  But to me (BW) this doesn't seem enough.  I think we
        % need to create a new lens file, right?  Potentially, we need to
        % adjust the accommodation, too.
        switch lower(val)
            case {'navarro'}
                obj.modelName = val;
                obj.recipe.set('retina distance',16.32);
                obj.recipe.set('retina radius',12);
            case { 'legrand'}
                obj.modelName = val;
                obj.recipe.set('retina distance',16.6);
                obj.recipe.set('retina radius',13.4);
            case {'arizona'}
                obj.modelName = val;
                obj.recipe.set('retina distance',16.713);
                obj.recipe.set('retina radius',13.4);
            otherwise
                % User defined name
                obj.modelName = val;
                fprintf('Custom model. Set the retina distance and radius\n');
        end
        
    case {'recipe'}
        % This is the iset3d recipe.
        obj.recipe = val;
        
    case 'lensdensity'
        % Adjust the lens pigment density
        obj.lensDensity = val;
        
    case 'usepinhole'
        % thisEye.set('use pinhole',true/false)
        %
        % Tells PBRT to use pinhole and generate a scene rather than using
        % the optics and generate an optical image. When true, the lens
        % will be replaced by a pinhole in the write() operation.
        obj.usePinhole = val;
        
    case 'useoptics'
        % Somes it just reads better to use
        %   thisEye.set('use optics',true/false)
        % instead of
        %   thisEye.set('use pinhole',false);
        obj.usePinhole = ~val;
                
    case 'fov'
        % We have a PPT about the parameters that need to be adjusted to
        % set the FOV for a realisticEye mode. The PowerPoint
        % (EyeballGeometry.pptx) is in the ISETBio/wiki/images directory.
        %
        % Setting the field of view amounts to setting the 'retina
        % semidiam' parameter. We figure out what it should be set to
        % here.
        %
        % The key equations are this - but you should look at the PPT.
        %
        %   fov = atand(semidiam/lens2chord)*2
        %   tand(fov/2) = semidiam/lens2chord
        %   semidiam = tand(fov/2)*lens2chord
        lens2chord  = obj.get('lens 2 chord','mm');
        semidiam    = tand(val/2)*lens2chord;
        radius      = obj.get('eye radius','mm');
        if semidiam >= radius
            error('Semidiam %f must be smaller than eyeball radius %f ',semidiam,radius);
        end
        obj.set('retina semidiam',semidiam);
      
    otherwise
        % If it is not an eyeSet parameter, it is probably an iset3d recipe
        % parameter. Send it in and hope for the best.
        obj.recipe.set(param,val,varargin{:});
end

end
