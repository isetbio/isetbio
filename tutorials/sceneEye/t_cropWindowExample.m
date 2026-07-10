%% t_cropWindowExample.m
%
% This demonstrates how to use the crop window functionality. Using the
% crop window let's you only render a specific portion of the scene. This
% can save you a lot of computation time.
% 
% Depends on: iset3d, isetbio, Docker
%
% TL ISETBIO Team, 2017  

%% Initialize
if piCamBio
    fprintf('%s: requires ISETBIO, not ISETCam\n',mfilename); 
    return;
end
    
ieInit;

%% Load a scene

scene3d = sceneEye('chessSet');
               
scene3d.fov = 30; 
scene3d.resolution = 128; 
scene3d.numRays = 128;
scene3d.numCABands = 0;
scene3d.accommodation = 1/0.4; % Accommodate to the pawn.

scene3d.name = 'chessSet_full';
oi = scene3d.render();
ieAddObject(oi);
oiWindow;

%% Pick a portion of the scene to re-render at a higher quality

rgb_full = oiGet(oi,'rgb');
x_full= scene3d.angularSupport;

% You can use getrect, but here we'll preset it.
if(scene3d.resolution == 128)
    % This rectangle is defined from the 128 resolution image.
    r = [77 63 20 20];
else
    error(['Incorrect resolution. Use getrect() to get a new',...
    ' rectangle for your resolution chosen.']);
end

% Show the rectangle we will render
figure();
image(rgb_full);
rectangle('Position',r,'EdgeColor','r','LineWidth',2)
xlabel('pixels'); ylabel('pixels');
axis image;

% Convert the rectangle to a crop window. According to pbrt-v3, the
% cropwindow is defined as follows:
%   The subregion of the image to render. The four values specified should
%   be fractions in the range [0,1], and they represent x_min, x_max,
%   y_min, and y_max, respectively. These values are in normalized device
%   coordinates, with (0,0) in the upper-left corner of the image.
cropwindow_px = [r(1) r(1)+r(3) r(2) r(2)+r(4)];
cropwindow_norm = cropwindow_px./scene3d.resolution;

% Set the crop window.
% The recipe contains all the instructions used to render. It's a structure
% primarily used in the ISET3d code. Typically when using sceneEye, we
% don't need to manipulate any values in the recipe. This is a special case
% where we have to.
scene3d.recipe.set('cropwindow',cropwindow_norm);

% Re-render
scene3d.name = 'chessSet_cropped_lowRes';
oi = scene3d.render();
ieAddObject(oi);
oiWindow;

%% Increase resolution
% The crop window has pretty low resolution. We can increase the resolution
% of scene3d to get a higher resolution cropped image.

% Let's assume we want a cropped window with resolution of 128 by 128.
desiredRes = [128 128];
currRes = [r(3) r(4)];

% Let's increase the resolution of the overall, full image so the cropped
% window will have our desired resolution.
scalingFactor = desiredRes./currRes;

% The full image has to be square, so we can only take one scaling value
scalingFactor = max(scalingFactor); 

scene3d.resolution = round(scene3d.resolution.*scalingFactor);

% Re-render
scene3d.name = 'chessSet_cropped_highRes';
oi = scene3d.render();
ieAddObject(oi);
oiWindow;

%% A warning on angular support

% sceneEye has a parameter for the angular support; this let's us plot the
% image in units of visual angle:
figure();
image(x_full,x_full,rgb_full); % These variables were saved earlier in the script.
axis image; xlabel('deg'); ylabel('deg');

% However, when using the crop window the angular support still matches the
% full image before being cropped, as evidenced by it's size. (Remember, we
% scaled the full image so that the crop window would have a higher
% resolution).
size(scene3d.angularSupport)

% Therefore you will have to calculate the the correct "cropped" angular
% support if you want to plot using units of visual angle.

%% END