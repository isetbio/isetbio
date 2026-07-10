%% t_vergenceAccomm.m
% 
% We can do a simple simulation of vergence and accommodation by having the
% eyes converge and accommodate to a "moving" red sphere.
%
% Depends on: iset3d, isetbio, Docker, RemoteDataToolbox
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end
ieInit;

%% Load scene
myScene = sceneEye('chessSet');

%% Set parameters

myScene.resolution = 128; 
myScene.numRays = 128;

ipd = 64e-3; % Average interpupillary distance

% We save the original parameters so we have them around even after moving
% and rotating the eye. 
originalPos = myScene.eyePos;
startingPos = originalPos + [0 0 0.1]; % This position looks a little better.

% Needed because we will objects to the 3d world. This text block contains
% the "original" scene.
originalWorld = myScene.recipe.world; 

%% Create binocular retinal images

% We'll move the red sphere to these distances along the z-axis (optical
% axis).
sphereZ = [-0.2 0.0 0.2]; % Along z-axis

for ii = 1:length(sphereZ)
    
    % [x y z]
    % We want the sphere to be at same height as the eye but at different
    % distances away from it.
    spherePos = [0 startingPos(2) sphereZ(ii)];
    
    % Move the two eyes apart
    leftEyePos = startingPos - [ipd/2 0 0];
    rightEyePos = startingPos + [ipd/2 0 0];
    
    % Make them point at the sphere
    myScene.eyeTo = spherePos;
    
    % Reset world (Important!)
    % We want to start from the original, vanilla version of the 3d scene
    % with no sphere.
    myScene.recipe.world = originalWorld;
    
    % Add the red sphere
    myScene.recipe = piAddSphere(myScene.recipe,...
        'rgb',[1 0 0],...
        'radius',0.005,...
        'location',spherePos);
    
    % Accommodate to the sphere
    dist = sqrt(sum((spherePos -leftEyePos).^2)); % in mm
    myScene.accommodation = 1/dist;
    
    % Render the left eye
    myScene.eyePos = leftEyePos;
    myScene.name = sprintf('leftEye_%0.2fm',sphereZ(ii));
    [oi, result] = myScene.render;
    vcAddAndSelectObject(oi);
    oiWindow;
    
    % Render the right eye
    myScene.eyePos = rightEyePos;
    myScene.name = sprintf('rightEye%0.2fm',sphereZ(ii));
    oi = myScene.render;
    vcAddAndSelectObject(oi);
    oiWindow;
    
end