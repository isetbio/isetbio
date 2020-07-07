function obj = eyeSet(obj,param,val,varargin)
% Set methods for dependent variables
%
% Synopsis
%
% Brief description
%
% Input
%
% Optional key/val parmeters
%
% Output
%
% Description
%
% See also
%  sceneEye

%%  Force param to lower case, no spaces\
param = ieParamFormat(param);

%%  Main switch statment
switch param
    case 'name'
        obj.name = val;
    case 'resolution'
        % obj.set('resolution',[x,y]);
        % Resolution means the number of spatial samples on the film plane.
        
        if numel(val) == 1, val = [val,val]; end
        obj.recipe.set('film resolution',val);

    case 'raysperpixel'
        % obj.set('rays per pixel',val);
        %
        % Number of samples set out from each film pixel
        obj.recipe.set('rays per pixel',val);
    case{'maxdepth','bounces','nbounces'}
        obj.recipe.set('n bounces',val);
        
    case 'pupildiameter'
        obj.recipe.set('pupil diameter',val);
    % case 'fov'  
        
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
            case {'custom'}
                % Use Navarro as a default.
                % Do we need to be able to change this later?
                obj.modelName = val;
                if(isempty(obj.retinaDistance))
                obj.recipe.set('retina distance',[]);
                end
                if(isempty(obj.retinaRadius))
                obj.recipe.set('retina radius',[]);
                end
        end
        
        % When the user toggles into debugMode, that indicates the lens
        % will be replaced by a pinhole in the write() phase.
    case 'debugmode'
        obj.debugMode = val;
        
        %{
        % Too many things were happening here.  I do not think we need this
        % any more.
        if(val)
            obj.modelName = 'perspective';
            % The camera will be changed to perspective in write(), so we
            % do nothing here.
        elseif(~val && strcmp(obj.modelName, 'none'))
            % Not debug mode. Put the navarro eye back in if there's not
            % already a model.  We should add in the accommodation
            % distance.
            obj.modelName = 'navarro';
            warning('Need to set the accom distance coming out of debug mode');
            obj.recipe.camera = piCameraCreate('humaneye');
        end
        %}
        
        % I want to put in this warning, but again MATLAB doesn't really like
        % this!
    case 'accommodation'
        % obj.set('accommodation',diopters);
        %
        obj.recipe.set('focal distance',1/val);
        
        if(strcmp(obj.modelName, 'Gullstrand'))
            warning(['Gullstrand eye has no accommodation modeling.', ...
                'Setting accommodation will do nothing.']);
        end
    case {'focaldistance'}
        % obj.set('focal distance',Meters)
        %
        obj.recipe.set('focal distance',1/val);
    case 'fov'
        % We specify our hope for the horizontal field of view
        obj.fov = val;
        
    case 'lensfile'
        % On creation, the lensFile is left empty
        % There should be a better way to do this right? I don't think I'm
        % doing this right. Maybe we need a seperate set function like we
        % do with render recipes?
        if(~isempty(val))
            % Make sure it's a valid file
            [p, ~, e] = fileparts(val);
            if(isempty(p))
                error('Lens file needs to be a full file path.');
            elseif(~strcmp(e, '.dat'))
                error('Lens file needs to be a .dat file.');
            end
            
            obj.modelName = 'custom';
            obj.recipe.set('lensFile',val);
        else
            obj.recipe.set('lensFile','');
        end
    case 'diffraction'
        % obj.set('diffraction',true/false);
        %
        obj.recipe.set('diffraction',val);
        
        % Eye position
    case 'lookat'
        % obj.set('look at',valStruct);
        %
        % Includes the from, to and up in a struct
        if isstruct(val) &&  isfield(val,'from') && isfield(val,'to')
            obj.recipe.lookAt = val;
        end
    case 'from'
        % obj.set('from',val);
        obj.recipe.lookAt.from = val;
    case 'to'
        obj.recipe.lookAt.to = val;
    case 'up'
        obj.recipe.lookAt.up = val;
end

end
