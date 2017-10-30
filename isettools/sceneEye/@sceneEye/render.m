function [oi, terminalOutput, outputFile] = render(obj, varargin)
%RENDER Render a scene3D object and return an optical image.
%   Given a scene3D object, we have all the information we need to
%   construct a PBRT file and render it. Therefore, this function does the
%   following:
%
%   1. Write out a new PBRT ([renderName].pbrt) in the working directory
%   2. Render using docker container
%   3. Load the output into an ISETBIO optical image, filling in the right
%   parameters with the scene information
%   4. Return the OI


%%  Given the scene3D object, we make adjustments to the renderecipeecipe 
recipe = obj.recipe;

% Apply any accommodation changes
recipe = setAccommodation(recipe,obj.accommodation,obj.workingDir);

% Film parameters
[recipe.film.name.value] = sprintf('%s.pbrt',obj.name);
[recipe.film.name.type] = 'string';
recipe.film.xresolution.value = obj.resolution;
recipe.film.yresolution.value = obj.resolution;

% Camera parameters
if(obj.debugMode)
    % Use a perspective camera with matching FOV instead of an eye.
    fov = struct('value',obj.fov,'type','float');
    recipe.camera = struct('type','Camera','subtype','perspective','fov',fov);
else
    recipe.camera.retinaDistance.value = obj.retinaDistance;
    recipe.camera.pupilDiameter.value = obj.pupilDiameter;
    recipe.camera.retinaDistance.value = obj.retinaDistance;
    recipe.camera.retinaRadius.value = obj.retinaRadius;
    recipe.camera.retinaSemiDiam.value = obj.retinaDistance*tand(obj.fov/2);
end

% Sampler
recipe.sampler.pixelsamples.value = obj.numRays;

% Integrator
recipe.integrator.maxdepth.value = obj.numBounces;

% Renderer
if(obj.numCABands == 0 || obj.numCABands == 1 || obj.debugMode)
    % No spectral rendering
    recipe.renderer = struct('type','Renderer','subtype','sampler');
else
    % Spectral rendering
    nWaveBands = struct('value',obj.numCABands,'type','integer');
    recipe.renderer = struct('type','Renderer', ...
        'subtype','spectralrenderer',...
        'nWaveBands',nWaveBands);
end

% Look At
if(isempty(obj.eyePos) || isempty(obj.eyeTo) || isempty(obj.eyeUp))
    error('Eye location missing!');
else
    recipe.lookAt = struct('from',obj.eyePos,'to',obj.eyeTo,'up',obj.eyeUp);
end

%% Write out the adjusted recipe into a PBRT file
pbrtFile = fullfile(obj.workingDir,strcat(obj.name,'.pbrt'));
recipe.outputFile = piWrite(recipe,pbrtFile,'overwrite',true);

%% Render the pbrt file using docker
[oi, outFile] = piRender(pbrtFile,'opticsType','lens');

%% Set OI parameters correctly:

% Scene distance. We set it to infinity, since it doesn't technically apply
% to the raytracing. 
oi = oiSet(oi, 'distance', Inf);

% This is the distance between the lens and the focal plane. This is not
% exactly the same as the focal length for the eye, but it's close. 
oi = oiSet(oi, 'optics focallength', obj.retinaDistance * 1e-3); 

oi = oiSet(oi,'optics fnumber',obj.retinaDistance/obj.pupilDiameter);
oi = oiSet(oi,'fov',obj.fov);
     
% Clear most of the default optics
% TODO: Is doing this okay? What's a better way to do this? 
oi.optics = opticsSet(oi.optics,'model','raytrace');
oi.optics = opticsSet(oi.optics,'name','PBRT Navarro Eye');
oi.optics.OTF = [];
oi.optics.lens = [];
oi.optics.offaxis = '';
oi.optics.vignetting = [];


end

