%% t_eyeIntro.m
%
% This brief tutorial is a (mainly) textual, introduction to modeling the
% optics of the human eye using ray-tracing in ISETBio and ISET3d.  
%
% The words and brief code in here describe the basic ideas.  Other
% tutorials (t_eye*) for more examples of parameters (with fewer words) to
% illustrate the ISETBio sceneEye usage.  Look through the ISET3d tutorials
% (t_piIntro*) to see more about general modeling with many other types of
% optical models.
% 
% To begin, you must have the Github repo ISET3d on your MATLAB path:
%
%   https://github.com/ISET/iset3d 
%
% as well as the Github repo isetbio on your path:
%
%   https://github.com/isetbio/isetbio
% 
% You must have Docker installed and running on your machine. You can find
% general instructions on docker here: 
%
%   https://www.docker.com/
%
% You can find the source code for pbrt-v3-spectral here:
% 
%   https://github.com/scienstanford/pbrt-v3-spectral
%
% Depends on: 
%    ISET3d, ISETBio, Docker
%
% See also:
%   YouTube videos:
%
  
%% The basic idea
%
% ISETBIO includes a set of tools for calculating the retinal image (also
% called the spectral irradiance at the retina). 
%
% Many of the tools were developed for relatively simple scenes, like the
% image on the central few degrees on a computer screen.  For such scenes
% the spectral point spread functions are enough for the central fovea.  
%
% If the image spans a larger field of view, then the spectral point spread
% functions across different field heights are necessary.
%
% For the much larger world of 3D natural images, we need additional tools.
% In part this is because distance from the eye's accommodative plane has
% an impact on the blur, and in part this is because 3D occlusion has an
% additional impact.
%
% To assist with calculating for 3D natural image models we use ray tracing
% methods. Specifically, we modified PBRT (Physically Based Ray Tracer) to
% enable us to calculate how a 3D scene would become a retinal image
% given that we had a model of the physiological optics.  PBRT is a widely
% used and admired ray tracer that has been validated, it is open-source,
% and it is very well-documented.
%
% PBRT is not Matlab, and thus we needed to develop a method for using it
% smoothly with ISETBio (and ISETCam).
%
% What we have done is place the PBRT code inside of a Docker container
% that can be called from within Matlab.  We have written a large number of
% Matlab tools that permit us to set up the PBRT scenes and the parameters
% of the physiological optics that transform the 3D scene into the retinal
% image.
%
% With this approach, you do not need to compile or install the source code
% in order to render images. You must only have docker installed and
% running on your computer. 
%
% The next few lines show the basic philosoph we take to running the 3D ray
% tracing code. There are many more tutorial scripts that go into detail,
% and we hope you find those useful.  This one just gets you going.
% Slightly.
%

%% Initialize ISETBIO

% We switch between ISETBio and ISETCam often. We use this little header to
% make sure the user is on the ISETBio path needed for sceneEye.
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end

% General ISET initialization
ieInit;

% Checks that you have docker installed and configured.
if ~piDockerExists, piDockerConfig; end

%% Render a fast, low quality retinal image

% We have several scenes that have been created specifically for ISETBio
% and ISET3d. 

% You can select and read the data from a PBRT scene by initializing a
% sceneEye object with the directory name of the PBRT files.  A collection
% of files is part of the ISET3d repository (iset3d/data/V3).
thisSE = sceneEye('letters at depth');

% theScene is a sceneEye object that has been initialized by reading the
% PBRT files in the relevant directory.
%
% To see all of PBRT files available for loading use the @recipe 'list'
% method:
thisSE.recipe.list;

% The first level of the sceneEye object class is quite simple, containing
% only a few parameters that are special to its function.
disp(thisSE)

% Almost all of the rendering parameters and the complexity of the scene
% are specified within the 'recipe'. That object (@recipe in ISET3d) is not
% simple.
%
% The @recipe class is used by ISET3d for all scenes, whether they are
% physiological optics models or not. The rendering software can handle
% many types of lenses and much more complex optical models, say with
% microlens arrays.  This is the top level of that class.  Notice that we
% retrieve the recipe from theScene using the 'get' method.
disp(thisSE.get('recipe'))

%%  Seeing the scene.

% It is often useful to have a quick look at the scene you have loaded. Ray
% tracing through a pinhole is much faster than through the optics. We take
% a quick look by setting the 'use pinhole' to true and then rendering.
% Because there are no optics, the result of the rendering is an ISETBio
% scene, rather than a spectral irradiance at the retina.

% Notice the use of 'set'
thisSE.set('use pinhole',true);

% We can tell PBRT the field of view of the scene.
thisSE.set('fov',40);

% Render it this way.   By default, the depth map is calculated, too.
% Notice that the rendering calls docker (twice).  The radiance and depth
% maps are calculated and read into the ISETBio scene structure.
scene = thisSE.render;

% Show the scene
sceneWindow(scene);

%% Rendering the PBRT scene

% Now tell PBRT to use the lens
thisSE.set('use optics',true);

% Set the units in mm (?) - default is true in PBRT realisticEye
thisSE.set('mmUnits', true);

% Set the scene to focus on the numbers at 200 mm
thisSE.set('accommodation',1/0.2);   % Diopters

% The lens rendering benefits from adding a more rays per pixel
thisSE.set('rays per pixel',192);

% Summarizing is always nice
thisSE.summary;

% Render and show.  We use 'oi' to refer to optical image, the spectral
% irradiance at the retina (also retinal irradiance).
oi = thisSE.render;
oiWindow(oi);

% The number at 200 is in good focus, while the others are not.

%% Set the accomodation to the 100 mm target

% Set the focal distance to 100 mm (0.1 m).  We will set 'accommodation'
% again, but you could set the 'focal distance' parameter if you prefer.
thisSE.set('accommodation',1/0.1);   % Diopters

% Summarize
thisSE.summary;

% Now we 
oi = thisSE.render;
oiWindow(oi);

% Now the number at 100 mm is in good focus, not 200 or 300.

%% Adjust the lens pigment density

% Other parameters, such as lens density are managed in the same way.  For
% example you can change the lens density to 0 this way.
thisSE.set('lens density',0);

thisSE.summary;

oi = thisSE.render;

oiWindow(oi);

% The image is no longer yellow.

%% More on the parameters

% There are many parameters of the rendering and eye model that can be
% changed. We always address this parameters through the 'get' and 'set'
% methods of the sceneEye class.  This is important because if we ever have
% to make a change to the parameters or calculations, that change can be
% implemented in the 'set' or 'get' function, rather than having to track
% down all the places that made explicit reference to the parameter.  All
% programmers know this - but we have observed that the idea is not alway
% implemented in practice.
%

% The parameters include values that control the rendering, control the
% eye, control objects in the scene.  At this point nearly all of the work
% we are doing on ISET3d - and it is a work in progress - involves better
% understanding and controlling the parameters needed to do scientific and
% engineering work.  We are not done - but we are better now than we were a
% year ago, so there is hope.

% We suggest that a next tutorial to try might be t_eyeFocalDistance, which
% illustrates chromatic aberration and focal distance for the Navarro
% physiological optics model.

%% END


