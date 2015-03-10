% s_opticsDepthScene
%
%  Illustrates an imperfect method for computing
%  and various rendering options.
%
%
% Copyright ImagEval Consultants, LLC, 2011

%% Initialize ISET Session
ieInit

%% Make a script/function to load the scene

% This file contains a scene with a depth map.  The scene was created
% outside of ISET as part of a project being done with Andy Lin.
load(fullfile(isetRootPath, 'data', 'scenes','piano3d.mat'));
scene = sceneSet(scene,'fov',3);

%%  Make optics with a little bigger pupil for depth of field effects

% Parameters for a large aperture
fNumber = 4;
pupilFactor = 3;  % When set to 1, this becomes diffraction limited.

oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics, 'otfmethod', 'custom');
optics = opticsSet(optics, 'model', 'ShiftInvariant');
optics = opticsSet(optics,'fnumber',fNumber);

% We only set focal length, not pupil diameter.  This trick keeps the f/#
% and adjusts the focal length so that the pupil becomes large
f = opticsGet(optics,'focal length');
optics = opticsSet(optics,'focal length',pupilFactor*f);

% Attach the optics to the oi and move on.
oi = oiSet(oi,'optics',optics);

%% These are the object distances for different defocus levels.

% We track depth edges over this defocus range
defocus = linspace(-1.2,0,7);
inFocusDepth = [1.5,100];
for ii=1:length(inFocusDepth)
    
    % This is the depth we would like to be in focus (m)
    thisFocusDepth = inFocusDepth(ii);
    
    % Find the depth edges so that the range of defocus is as above, though
    % it is centered around a depth of inFocusDepth.  To achieve this the
    % imageDist will not be in the focal plane.
    [depthEdges, imageDist, oDefocus] = oiDepthEdges(oi,defocus,thisFocusDepth);
    
    % Get the scene depth map, blur it, reattach it to the scene
    oMap  = sceneGet(scene,'depth map');
    sceneDepthRange = [depthEdges(1),10];
    oMap  = ieScale(oMap,sceneDepthRange(1),sceneDepthRange(2));
    blurSize = 2; supportSize = [5 5];
    g = fspecial('gaussian',supportSize,blurSize); 
    oMap = conv2(oMap,g,'same');
    scene = sceneSet(scene,'depth map',oMap);
    % vcAddAndSelectObject(scene); sceneWindow;
    
    % Main function.
    cAberration = [];
    displayFlag = 0;
    [oiD, D] = oiDepthCompute(oi,scene,imageDist,depthEdges,cAberration,displayFlag);
    % for ii=1:length(oiD), vcAddAndSelectObject(oiD{ii}); end; oiWindow
    
    % Combine them and show them in the window.
    oi = oiDepthCombine(oiD,scene,depthEdges);
    pupil = opticsGet(optics,'pupil radius','mm');
    oi = oiSet(oi,'name',sprintf('Focus-%.1fm',thisFocusDepth));
    
    vcAddAndSelectObject(oi);  oiWindow
    
end

%%
imageMultiview('oi',[1 2], true);

%% End