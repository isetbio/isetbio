%% Figure generation for RealityLabs 2021 application
%
% ISETBio-ISET3D

thisDir = '/Volumes/GoogleDrive/My Drive/Grants/2018 Occulus (Hillis)'; 
chdir(thisDir)

%% The cone mosaic at different eccentricities 

% A long strip along the horizontal axis
cm = cMosaic('positionDegs',[5 0],'sizeDegs',[10 2]);
cm.visualize;
title('');

%% Rendering of a uniform field


% We can see the 'breaks' in the mosaic in this image
scene = sceneCreate('uniform ee');
scene = sceneInterpolate(scene,10);
scene = sceneSet(scene,'hfov',10);
sceneWindow(scene);

oi = oiCreate;
oi = oiCompute(oi,scene);
[allE, allNoisyE] = cm.compute(oi);
cm.plot('excitations horizontal line',allE,'ydeg',0);
cm.plot('excitations',allE);

%% Rendering of a grating

% The breaks in the mosaic are less obvious and we can make it rectangular
% more easily.
hp = harmonicP; hp.freq = 5; hp.row = 128; hp.col = 512;
scene = sceneCreate('harmonic',hp);
scene = sceneSet(scene,'hfov',10);
sceneWindow(scene);

oi = oiCreate;
oi = oiCompute(oi,scene);
[allE, allNoisyE] = cm.compute(oi);
uData = cm.plot('excitations horizontal line',allNoisyE,'ydeg',0);
cm.plot('excitations',allE);
cm.visualize;
title('');

% Excitations per degree is nearly constant, but that is probably only
% because we omitted the rods.  The number would go down in reality.
excitations = squeeze(uData.roiE);
for ii=0:2:8
    l = (uData.pos(:,1) > ii) & (uData.pos(:,1) < ii + 2);
    sum(excitations(l))
end

%% Rendering of letters
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end
ieInit;
if ~piDockerExists, piDockerConfig; end

%% Here are the World positions of the letters in the scene

% The units are in meters
toA = [-0.0486     0.0100     0.5556];
toB = [  0         0.0100     0.8333];
toC = [ 0.1458     0.0100     1.6667];

%% Show the scene

% This is rendered using a pinhole so the rendering is fast.  It has
% infinite depth of field (no focal distance).
thisSE = sceneEye('letters at depth','eye model','legrand');
% thisSE.summary;

% Position the eye off to the side so we can see the 3D easily
from = [0.25,0.3,-1.3];
thisSE.set('from',from);

% Look at the position with the 'B'.  The values for each of the letters
% are included above.
thisSE.set('to',toB);
thisSE.set('from',from - [0.5 0 0.7]);

% Reduce the rendering noise by using more rays. 
thisSE.set('rays per pixel',32);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('spatial samples',512);  

% Have a quick check with the pinhole
thisSE.set('use pinhole',true);

% thisSE.get('object distance')   % Default is 2.1674
% If we make it further, we can narrow the FOV, I think
% thisSE.set('object distance',6);
% thisSE.set('fov',6);

% Given the distance from the scene, this FOV captures everything we want
thisSE.set('fov',6);             % Degrees

letterA = piAssetSearch(thisSE.recipe,'object name','_A');
letterB = piAssetSearch(thisSE.recipe,'object name','_B');
letterC = piAssetSearch(thisSE.recipe,'object name','_C');

%{
letterA = '001_A_O';
letterB = '001_B_O';
letterC = '001_C_O';
%}

bPos = thisSE.get('asset',letterB,'world position');
cPos = thisSE.get('asset',letterC,'world position');

% The fraction of the distance from C to B.
thisSE.set('asset',letterC,'world position', cPos + 0.5*(bPos - cPos));

%% Insert the materials we want
piMaterialsInsert(thisSE.recipe,'groups',{'glass','diffuse','wood','brick'});

% For glass and metal, multiple bounces is needed.
thisSE.set('nbounces',5);

thisSE.get('print materials')

%%
thisSE.set('asset', letterA, 'material name', 'White');
thisSE.set('asset', letterB, 'material name', 'glass');
thisSE.set('asset', letterC, 'material name', 'diffuse-gray');

% Now make an image texture out of grass
% Scale is cumulative.  So ... don't run this multiple times
ground = piAssetSearch(thisSE.recipe,'object name','Ground');
thisSE.set('asset',ground,'scale',[0.15 .3 1]);
thisSE.set('asset', ground, 'material name', 'wood-light-large-grain');

% Do the wall
wall = piAssetSearch(thisSE.recipe,'object name','Wall');
thisSE.set('asset',wall,'scale',[0.15 1 0.15]);
thisSE.set('asset', wall, 'material name', 'brickwall001');
thisSE.recipe.show('objects');

%% Render the scene at different camera positions.  
% 
% Add those in to make some point about stereo.
from = thisSE.get('from');
thisSE.recipe.set('render type',{'radiance'}); 

%%
thisSE.set('from',from - [0.06 0.0 0]);
scene = thisSE.render; 
sceneWindow(scene);

%%
thisSE.set('from',from + [0.06 0.0 0]);
scene = thisSE.render; 
sceneWindow(scene);

%% Use the eye model.  Requires CPU, however.

% Turn off the pinhole.  The model eye (by default) is the Navarro model.
thisSE.set('use pinhole',false);

% We turn on chromatic aberration.  That slows down the calculation, but
% makes it more accurate and interesting.  We oftens use only 8 spectral
% bands for speed and to get a rought sense. You can use up to 31.  It is
% slow, but that's what we do here because we are only rendering once. When
% the GPU work is completed, this will be fast!
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

% Find the distance to the object
oDist = thisSE.get('object distance');

% This is the distance to the B and we set our accommodation to that.
thisSE.set('focal distance',oDist);  

% Reduce the rendering noise by using more rays. 
thisSE.set('rays per pixel',256);      

% Increase the spatial resolution by adding more spatial samples.
thisSE.set('spatial samples',512);     

thisSE.set('from',from);
oi = thisSE.render('render type','radiance');
oiWindow(oi);

% This takes longer than the pinhole rendering, so we do not bother with
% the depth.
thisSE.set('from',from - [0.06,0,0]);
oi = thisSE.render('render type','radiance');
oiWindow(oi);

thisSE.set('from',from + [0.06,0,0]);
oi = thisSE.render('render type','radiance');
oiWindow(oi);

save(fullfile(thisDir,'oiLetters'),'oi');

%% Render the 3D letters on the mosaic


% load('cm1x1.mat','cm');
load('cm5x5.mat','cm');

% A long strip along the horizontal axis
% cm = cMosaic('positionDegs',[0 0],'sizeDegs',[1 1]);
cm.visualize;
title('');

oi = vcGetObject('oi');

% save('cm5x5','cm');
[allE, noisyE] = cm.compute(oi);
cm.plot('excitations',allE);
title('');

%% Rendering of stereo pairs

