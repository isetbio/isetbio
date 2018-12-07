function [objNew] = write(obj, varargin)
% Write out the sceneEye object into a pbrt file that will later be
% rendered. This function is a substep of sceneEye.render. Typically a user
% will not run this directly, but rather it will be run within the render
% function.
%
% Syntax:
%   [success] = write(obj, [varargin])
%
% Description:
%	 Given a sceneEye object, we have all the information we need to
%    construct a PBRT file and render it. Therefore, this function reads
%    and interprets the parameters given in the sceneEye object and writes
%    them out into an PBRT file. This file will later be rendered in
%    sceneEye.render.
%
% Inputs:
%    obj            - The scene3D object to render
%    varargin       - (Optional) Other key/value pair arguments
%
% Outputs:
%   objNew          - the object may have been modified in the processing
%                     below. We return this modified version.
%

%% Make a copy of the current object
% We will render this copy, since we may make changes to certain parameters
% before rendering (i.e. in th eccentricity calculations) but we don't want
% these changes to show up original object given by the user.
objNew = copy(obj);
objNew.recipe = copy(obj.recipe);

%% Make some eccentricity calculations

% To render an image centered at a certain eccentricity without having
% change PBRT, we do the following:
% 1. Change the film size and resolution so that renders a larger image
% that encompasses the desired eccentricity (tempWidth/tempHeight)
% 2. Insert a "crop window" PBRT parameter to only render the window
% centered at the desired eccentricity with the desired film
% diagonal/resolution.

ecc = objNew.eccentricity;

% I was having many bugs with my eccentricity code, so for now I've removed
% it for now. Ideally we do all the right calculations shown above and then
% use scene3d.recipe.set('cropwindow',[x1 x2 y1 y2]); and then carefully
% reset the angular support as well...
if(ecc ~= [0 0])
    warning('Eccentricity is currently not implemented. Setting to zero.')
    ecc = [0 0];
end

%% Given the sceneEye object, we make all other adjustments needed to the recipe
recipe = objNew.recipe;

% Depending on the eye model, set the lens file appropriately
switch objNew.modelName
    case {'Navarro','navarro'}
        % Apply any accommodation changes
        if(isempty(objNew.accommodation))
            objNew.accommodation = 5;
            warning('No accommodation! Setting to 5 diopters.');
        end
        
        % This function also writes out the Navarro lens file
        recipe = setNavarroAccommodation(recipe, objNew.accommodation,...
                                         objNew.workingDir);
        
    case {'Gullstrand','gullstrand'}
        
        % Gullstrand eye does not have accommodation (not yet at least), so
        % for now all we need to do is write out the lens file.
        
        lensFile = 'gullstrand.dat';
        writeGullstrandLensFile(fullfile(objNew.workingDir, lensFile));
        fprintf('Wrote out a new lens file: \n')
        fprintf('%s \n \n', fullfile(objNew.workingDir, lensFile));
        
        objNew.recipe.camera.lensfile.value = fullfile(objNew.workingDir, lensFile);
        objNew.recipe.camera.lensfile.type = 'string';   
    
    case{'Arizona','arizona'}
        
        if(isempty(objNew.accommodation))
            objNew.accommodation = 5;
            warning('No accommodation! Setting to 5 diopters.');
        end
        
        % This function also writes out the Arizona lens file.
        recipe = setArizonaAccommodation(recipe, objNew.accommodation,...
                                         objNew.workingDir);
                                     
end


% Film parameters
recipe.film.xresolution.value = objNew.resolution;
recipe.film.yresolution.value = objNew.resolution;

% Camera parameters
if(objNew.debugMode)
    % Use a perspective camera with matching FOV instead of an eye.
    fov = struct('value', objNew.fov, 'type', 'float');
    recipe.camera = struct('type', 'Camera', 'subtype', 'perspective', ...
        'fov', fov);
    if(objNew.accommodation ~= 0)
        warning(['Setting perspective camera focal distance to %0.2f dpt '...
            'and lens radius to %0.2f mm'],...
            objNew.accommodation,objNew.pupilDiameter);
        recipe.camera.focaldistance.value = 1/objNew.accommodation;
        recipe.camera.focaldistance.type = 'float';
        
        recipe.camera.lensradius.value = (objNew.pupilDiameter/2)*10^-3;
        recipe.camera.lensradius.type = 'float';
    end
else
    recipe.camera.retinaDistance.value = objNew.retinaDistance;
    recipe.camera.pupilDiameter.value = objNew.pupilDiameter;
    recipe.camera.retinaDistance.value = objNew.retinaDistance;
    recipe.camera.retinaRadius.value = objNew.retinaRadius;
    recipe.camera.retinaSemiDiam.value = objNew.retinaDistance ...
        * tand(objNew.fov / 2);
    if(strcmp(objNew.sceneUnits,'m'))
        recipe.camera.mmUnits.value = 'false';
        recipe.camera.mmUnits.type = 'bool';
    end
    if(objNew.diffractionEnabled)
        recipe.camera.diffractionEnabled.value = 'true';
        recipe.camera.diffractionEnabled.type = 'bool';
    end
end

% Sampler
recipe.sampler.pixelsamples.value = objNew.numRays;

% Integrator
recipe.integrator.maxdepth.value = objNew.numBounces;
recipe.integrator.maxdepth.type = 'integer';

% Renderer
if(objNew.numCABands == 0 || objNew.numCABands == 1 || objNew.debugMode)
    % No spectral rendering
    recipe.integrator.subtype = 'path';
else
    % Spectral rendering
    numCABands = struct('value', objNew.numCABands, 'type', 'integer');
    recipe.integrator = struct('type', 'Integrator', ...
        'subtype', 'spectralpath', ...
        'numCABands', numCABands);
end

% Look At
if(isempty(objNew.eyePos) || isempty(objNew.eyeTo) || isempty(objNew.eyeUp))
    error('Eye location missing!');
else
    recipe.lookAt = struct('from', objNew.eyePos, 'to', objNew.eyeTo, ...
        'up', objNew.eyeUp);
end

% If there was a crop window, we have to update the angular support that
% comes with sceneEye
% We can't do this right now because the angular support is a dependent
% variable. How to overcome this?
%{
currAngSupport = obj.angularSupport;
cropWindow = recipe.get('cropwindow');
cropWindowR = cropWindow.*obj.resolution;
cropWindowR = [cropWindowR(1) cropWindowR(3) ...
    cropWindowR(2)-cropWindowR(1) cropWindowR(4)-cropWindowR(3)];
[X,Y] = meshgrid(currAngSupport,currAngSupport);
X = imcrop(X,cropWindowR);
Y = imcrop(Y,cropWindowR);
% Assume square optical image for now, but we should probably change
% angularSupport to have both x and y direction.
objNew.angularSupport = X(1,:); 
%}

%% Write out the adjusted recipe into a PBRT file
pbrtFile = fullfile(objNew.workingDir, strcat(objNew.name, '.pbrt'));
recipe.set('outputFile',pbrtFile);
if(strcmp(recipe.exporter,'C4D'))
    piWrite(recipe, 'overwritepbrtfile', true, 'overwritelensfile', false, ...
        'overwriteresources', false,'creatematerials',true);
else
    piWrite(recipe, 'overwritepbrtfile', true, 'overwritelensfile', false, ...
        'overwriteresources', false);
end
obj.recipe = recipe; % Update the recipe.

end
