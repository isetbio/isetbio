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
% Notes:
%    * TODO: Is clearing most of the default optics okay? Determine if
%      there is a better way to accomplish this.
%    * TODO: Is it possible to remove the need to specify the type and
%      value for the chromaticAberrationEnabled?
%

%%  Given the scene3D object, we make adjustments to the renderecipeecipe 
recipe = obj.recipe;

% Apply any accommodation changes
if(isempty(obj.accommodation))
    obj.accommodation = 5;
    warning('No accommodation! Setting to 5 diopters.');
end
recipe = setAccommodation(recipe, obj.accommodation, obj.workingDir);

% Film parameters
recipe.film.xresolution.value = obj.resolution;
recipe.film.yresolution.value = obj.resolution;

% Camera parameters
if(obj.debugMode)
    % Use a perspective camera with matching FOV instead of an eye.
    fov = struct('value', obj.fov, 'type', 'float');
    recipe.camera = struct('type', 'Camera', 'subtype', 'perspective', ...
        'fov', fov);
else
    recipe.camera.retinaDistance.value = obj.retinaDistance;
    recipe.camera.pupilDiameter.value = obj.pupilDiameter;
    recipe.camera.retinaDistance.value = obj.retinaDistance;
    recipe.camera.retinaRadius.value = obj.retinaRadius;
    recipe.camera.retinaSemiDiam.value = obj.retinaDistance ...
        * tand(obj.fov / 2);
end

% Sampler
recipe.sampler.pixelsamples.value = obj.numRays;

% Integrator
recipe.integrator.maxdepth.value = obj.numBounces;

% Renderer
if(obj.numCABands == 0 || obj.numCABands == 1 || obj.debugMode)
    % No spectral rendering
    recipe.renderer = struct('type', 'Renderer', 'subtype', 'sampler');
else
    % Spectral rendering
    % TODO: Can we get rid of needing this?
    recipe.camera.chromaticAberrationEnabled.value = 'true';
    recipe.camera.chromaticAberrationEnabled.type = 'bool';
    
    nWaveBands = struct('value', obj.numCABands, 'type', 'integer');
    recipe.renderer = struct('type', 'Renderer', ...
        'subtype', 'spectralrenderer', ...
        'nWaveBands', nWaveBands);
end

% Look At
if(isempty(obj.eyePos) || isempty(obj.eyeTo) || isempty(obj.eyeUp))
    error('Eye location missing!');
else
    recipe.lookAt = struct('from', obj.eyePos, 'to', obj.eyeTo, ...
        'up', obj.eyeUp);
end

%% Write out the adjusted recipe into a PBRT file
pbrtFile = fullfile(obj.workingDir, strcat(obj.name, '.pbrt'));
recipe.outputFile = pbrtFile;
piWrite(recipe, 'overwritepbrtfile', true, 'overwriteresources', false, ...
    'overwritelensfile', false);

%% Render the pbrt file using docker
[ieObject, terminalOutput] = piRender(recipe);

%% Set OI parameters correctly:
if(~obj.debugMode)
    % Scene distance. We set it to infinity, since it doesn't technically
    % apply to the raytracing.
    ieObject = oiSet(ieObject, 'distance', Inf);
    
    % This is the distance between the lens and the focal plane. This is
    % not exactly the same as the focal length for the eye, but it's close.
    ieObject = oiSet(ieObject, 'optics focallength', ...
        obj.retinaDistance * 1e-3);
    ieObject = oiSet(ieObject, 'optics fnumber', ...
        obj.retinaDistance / obj.pupilDiameter);
    ieObject = oiSet(ieObject, 'fov', obj.fov);
    
    % Clear most of the default optics
    % TODO: Is doing this okay? What's a better way to do this?
    ieObject.optics = opticsSet(ieObject.optics, 'model', 'raytrace');
    ieObject.optics = opticsSet(ieObject.optics, 'name', ...
        'PBRT Navarro Eye');
    ieObject.optics.OTF = [];
    ieObject.optics.lens = [];
    ieObject.optics.offaxis = '';
    ieObject.optics.vignetting = [];
end

end
