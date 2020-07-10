function obj = eyeSet(obj,param,val,varargin)
% Set methods for dependent variables
%
% Synopsis
%   obj = eyeSet(obj,param,val,varargin)
%
% Brief description
%  Set the iset3d recipe and sceneEye parameters through this call.
%
% Input
%  obj:   sceneEye class
%  param:
%  val:
% 
% Optional key/val parmeters
%  N/A
%
% Output
%   obj:  Modified sceneEye class
%
% Description
%
%
% See also
%  sceneEye

%%  Force param to lower case, no spaces\
param = ieParamFormat(param);

%%  Main switch statment
switch param
    case 'name'
        obj.name = val;

    case 'modelname'
        % When we set the eye model, we need to change the retina distance and
        % radius.  
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
        
        % When the user toggles into debugMode, that indicates the lens
        % will be replaced by a pinhole in the write() phase.
    case 'usepinhole'
        obj.usePinhole = val;
                
    case 'fov'
        % We specify our hope for the horizontal field of view.  We have a
        % PPT about the various parameters needed here.  The PPTX is in the
        % wiki/images directory.
        obj.fov = val;
      
    otherwise
        % Rather than a sceneEye set, this is probably an iset3d recipe
        % set.  So send it in and hoe for the best.
        obj.recipe.set(param,val,varargin{:});
end

end
