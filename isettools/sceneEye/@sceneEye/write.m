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

%% Make some eccentricity calculations

% To render an image centered at a certain eccentricity without having
% change PBRT, we do the following:
% 1. Change the film size and resolution so that renders a larger image
% that encompasses the desired eccentricity (tempWidth/tempHeight)
% 2. Insert a "crop window" PBRT parameter to only render the window
% centered at the desired eccentricity with the desired film
% diagonal/resolution.

ecc = objNew.eccentricity;

% This section of the code has not been thoroughly finished/debugged, so
% let's put out a warning.
if(ecc ~= [0 0])
    warning('Eccentricity calculations have not been debugged. Use at your own risk!')
end

% Given a point at a certain eccentricitity [ecc(1) ecc(2)], what is
% the minimum FOV the rendered image needs to have in order to
% encompass the given point?
tempWidth = 2*obj.retinaDistance*tand(abs(ecc(1))) + obj.width;
tempHeight = 2*obj.retinaDistance*tand(abs(ecc(2))) + obj.height;
fovHoriz = 2*atand(tempWidth/(2*obj.retinaDistance));
fovVert = 2*atand(tempHeight/(2*obj.retinaDistance));
objNew.fov = max(fovHoriz,fovVert); 

% Center of image in mm, given desired ecc
centerX = obj.retinaDistance*tand(ecc(1));
centerY = obj.retinaDistance*tand(ecc(2));

% Boundaries of crop window in mm
% (Use original width and height!)
left = centerX - obj.width/2;
right = centerX + obj.width/2;
bottom = centerY + obj.height/2;
top = centerY - obj.height/2;

% Convert (0,0) to top left corner (normalized device coordinates) instead
% of center
tempSize = 2*objNew.retinaDistance*tand(objNew.fov/2); % Side length of large FOV
left_ndc = left + tempSize/2;
right_ndc = right + tempSize/2;
top_ndc = top + tempSize/2;
bottom_ndc = bottom + tempSize/2;
ndcWindow = [left_ndc right_ndc top_ndc bottom_ndc];

% Convert to ratio
cropWindowEcc = ndcWindow./tempSize;

% Since we'll be cropping the large image down to the desired
% eccentricity, we have to increase the rendered resolution.
tempResolution = objNew.resolution/(cropWindowEcc(2)-cropWindowEcc(1));
objNew.resolution = round(tempResolution);

% DEBUG
%{
    fprintf('*** DEBUG *** \n')
    fprintf('Original FOV: %0.2f \n',obj.fov);
    fprintf('New FOV: %0.2f \n',objNew.fov);
    fprintf('Original width: %0.2f \n', obj.width);
    fprintf('New width: %0.2f \n',objNew.width);
    fprintf('Original resolution: %0.2f \n',obj.resolution);
    fprintf('New resolution: %0.2f \n',objNew.resolution);
    fprintf('Crop window: [%0.2f %0.2f %0.2f %0.2f] \n',cropWindow);
    fprintf('*** DEBUG *** \n')
%}


%% Given the sceneEye object, we make all other adjustments needed to the recipe
recipe = objNew.recipe;

% Apply any accommodation changes
if(isempty(objNew.accommodation))
    objNew.accommodation = 5;
    warning('No accommodation! Setting to 5 diopters.');
end
recipe = setAccommodation(recipe, objNew.accommodation, objNew.workingDir);

% Film parameters
recipe.film.xresolution.value = objNew.resolution;
recipe.film.yresolution.value = objNew.resolution;

% Camera parameters
if(objNew.debugMode)
    % Use a perspective camera with matching FOV instead of an eye.
    fov = struct('value', objNew.fov, 'type', 'float');
    recipe.camera = struct('type', 'Camera', 'subtype', 'perspective', ...
        'fov', fov);
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

% Crop window

% Crop window and eccentricity can conflicting values. Let's resolve
% that (messily) here. We may rethink the eccentricity calculation in
% the future:
% If there is no eccentricity set, use whatever crop window was in the
% structure originally. If there is eccentricity, use the updated
% cropwindow. 
if(ecc == [0 0])
    % Do nothing
else
    recipe.film.cropwindow.value = cropWindowEcc;
    recipe.film.cropwindow.type = 'float';
end


%% Write out the adjusted recipe into a PBRT file
pbrtFile = fullfile(objNew.workingDir, strcat(objNew.name, '.pbrt'));
recipe.outputFile = pbrtFile;
piWrite(recipe, 'overwritepbrtfile', true, 'overwritelensfile', false, ...
    'overwriteresources', false);
obj.recipe = recipe; % Update the recipe.

end
