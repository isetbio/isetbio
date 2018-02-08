function [ieObject, terminalOutput] = render(obj, varargin)
% Render a scene3D object and return an optical image.
%
% Syntax:
%   [ieObject, terminalOutput, outputFile] = render(obj, [varargin])
%
% Description:
%	 Given a scene3D object, we have all the information we need to
%    construct a PBRT file and render it. Therefore, this function does the
%    following:
%       1. Write out a new PBRT ([renderName].pbrt) in working directory
%       2. Render using docker container
%       3. Load the output into an ISETBIO optical image, filling in the
%          right parameters with the scene information
%       4. Return the OI
%
% Inputs:
%    obj            - The scene3D object to render
%    varargin       - (Optional) Other key/value pair arguments
%
% Outputs:
%    ieObject       - The Optical Image object
%    terminalOutput - Terminal output
%

%% Make a copy of the current object
% We will render this copy, since we may make changes to certain parameters
% before rendering (i.e. in th eccentricity calculations) but we don't want
% these changes to show up original object given by the user.
objCopy = copy(obj);

%% Make some eccentricity calculations

% To render an image centered at a certain eccentricity without having
% change PBRT, we do the following:
% 1. Change the film size and resolution so that renders a larger image
% that encompasses the desired eccentricity (tempWidth/tempHeight)
% 2. Insert a "crop window" PBRT parameter to only render the window
% centered at the desired eccentricity with the desired film
% diagonal/resolution.

ecc = objCopy.eccentricity;

% Given a point at a certain eccentricitity [ecc(1) ecc(2)], what is
% the minimum FOV the rendered image needs to have in order to
% encompass the given point?
tempWidth = 2*obj.retinaDistance*tand(abs(ecc(1))) + obj.width;
tempHeight = 2*obj.retinaDistance*tand(abs(ecc(2))) + obj.height;
fovHoriz = 2*atand(tempWidth/(2*obj.retinaDistance));
fovVert = 2*atand(tempHeight/(2*obj.retinaDistance));
objCopy.fov = max(fovHoriz,fovVert); 

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
tempSize = 2*objCopy.retinaDistance*tand(objCopy.fov/2); % Side length of large FOV
left_ndc = left + tempSize/2;
right_ndc = right + tempSize/2;
top_ndc = top + tempSize/2;
bottom_ndc = bottom + tempSize/2;
ndcWindow = [left_ndc right_ndc top_ndc bottom_ndc];

% Convert to ratio
cropWindow = ndcWindow./tempSize;

% Since we'll be cropping the large image down to the desired
% eccentricity, we have to increase the rendered resolution.
tempResolution = objCopy.resolution/(cropWindow(2)-cropWindow(1));
objCopy.resolution = round(tempResolution);

% DEBUG
%{
    fprintf('*** DEBUG *** \n')
    fprintf('Original FOV: %0.2f \n',obj.fov);
    fprintf('New FOV: %0.2f \n',objCopy.fov);
    fprintf('Original width: %0.2f \n', obj.width);
    fprintf('New width: %0.2f \n',objCopy.width);
    fprintf('Original resolution: %0.2f \n',obj.resolution);
    fprintf('New resolution: %0.2f \n',objCopy.resolution);
    fprintf('Crop window: [%0.2f %0.2f %0.2f %0.2f] \n',cropWindow);
    fprintf('*** DEBUG *** \n')
%}


%% Given the sceneEye object, we make all other adjustments needed to the recipe
recipe = objCopy.recipe;

% Apply any accommodation changes
if(isempty(objCopy.accommodation))
    objCopy.accommodation = 5;
    warning('No accommodation! Setting to 5 diopters.');
end
recipe = setAccommodation(recipe, objCopy.accommodation, objCopy.workingDir);


% Film parameters
recipe.film.xresolution.value = objCopy.resolution;
recipe.film.yresolution.value = objCopy.resolution;

% Camera parameters
if(objCopy.debugMode)
    % Use a perspective camera with matching FOV instead of an eye.
    fov = struct('value', objCopy.fov, 'type', 'float');
    recipe.camera = struct('type', 'Camera', 'subtype', 'perspective', ...
        'fov', fov);
else
    recipe.camera.retinaDistance.value = objCopy.retinaDistance;
    recipe.camera.pupilDiameter.value = objCopy.pupilDiameter;
    recipe.camera.retinaDistance.value = objCopy.retinaDistance;
    recipe.camera.retinaRadius.value = objCopy.retinaRadius;
    recipe.camera.retinaSemiDiam.value = objCopy.retinaDistance ...
        * tand(objCopy.fov / 2);
end

% Sampler
recipe.sampler.pixelsamples.value = objCopy.numRays;

% Integrator
recipe.integrator.maxdepth.value = objCopy.numBounces;

% Renderer
if(objCopy.numCABands == 0 || objCopy.numCABands == 1 || objCopy.debugMode)
    % No spectral rendering
    recipe.renderer = struct('type', 'Renderer', 'subtype', 'sampler');
else
    % Spectral rendering
    nWaveBands = struct('value', objCopy.numCABands, 'type', 'integer');
    recipe.renderer = struct('type', 'Renderer', ...
        'subtype', 'spectralrenderer', ...
        'nWaveBands', nWaveBands);
end

% Look At
if(isempty(objCopy.eyePos) || isempty(objCopy.eyeTo) || isempty(objCopy.eyeUp))
    error('Eye location missing!');
else
    recipe.lookAt = struct('from', objCopy.eyePos, 'to', objCopy.eyeTo, ...
        'up', objCopy.eyeUp);
end

% Crop window
if(exist('cropWindow','var'))
    recipe.film.cropwindow.value = cropWindow;
    recipe.film.cropwindow.type = 'float';
end

%% Write out the adjusted recipe into a PBRT file
pbrtFile = fullfile(objCopy.workingDir, strcat(objCopy.name, '.pbrt'));
recipe.outputFile = pbrtFile;
piWrite(recipe, 'overwritepbrtfile', true, 'overwritelensfile', false, ...
    'overwriteresources', false);

%% Render the pbrt file using docker
[ieObject, terminalOutput] = piRender(recipe);

%% Set OI parameters correctly:
if(~obj.debugMode)
    % Scene distance. We set it to infinity, since it doesn't technically
    % apply to the raytracing.
    ieObject = oiSet(ieObject, 'distance', Inf);
    
    % For the optics focal length, we use the distance between the back of
    % the lens and the retina. Although it is not exactly the same as the
    % focal length for the eye (which also changes with accommodation), it
    % should be close enough for our purposes.
    ieObject = oiSet(ieObject, 'optics focallength', ...
        obj.retinaDistance * 1e-3);
    ieObject = oiSet(ieObject, 'optics fnumber', ...
        obj.retinaDistance / obj.pupilDiameter);
    ieObject = oiSet(ieObject, 'fov', obj.fov);
    
    % Clear default optics that do not apply to the ray-traced optical
    % image. We may want to add these in in the future.
    ieObject.optics = opticsSet(ieObject.optics, 'model', 'raytrace');
    ieObject.optics = opticsSet(ieObject.optics, 'name', ...
        'PBRT Navarro Eye');
    ieObject.optics.OTF = [];
    ieObject.optics.lens = [];
    ieObject.optics.offaxis = '';
    ieObject.optics.vignetting = [];
end

end
