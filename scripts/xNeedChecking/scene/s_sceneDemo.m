%% s_sceneDemo
%
%  Tutorial introduction to scripting, applied to a scene 
%
% The purpose of the script is to illustrate some of the simple scripting
% methods available in ISET.
%
% This script illustrates how to create a scene of the Macbeth ColorChecker illuminated with D65 
% and to display a downsampled (3 color channels) representation of the scene in the Display window.
% Various properties of the scene are read using
% sceneGet.  Some properties are set.  
%
% Then we create a frequency-orientation test target, extract the data, and
% make a plot of the luminance (cd/m2) across the bottom row of that
% target.
%
% see also s_sceneFromMultispectral, s_sceneFromRGB
%
% Copyright ImagEval Consultants, LLC, 2003.
%
%% sceneCreate
% To create a simple spectral scene
% of a Macbeth Chart  under a D65 illuminant, we use 
%
sceneMacbethD65 = sceneCreate('macbethd65');

%% sceneWindow
% To place the scene data in the window, add the scene object to the session and select it.
%
vcAddAndSelectObject('scene',sceneMacbethD65);
%
% Then bring up the scene window.  You can interact with the scene through this window
sceneWindow;
%
%% sceneGet
% To manipulate the data in a scene, you can extract variables
sceneGet(sceneMacbethD65,'meanLuminance')
%%
% Image of the luminance map of the Macbeth
luminance = sceneGet(sceneMacbethD65,'luminance');
vcNewGraphWin;
imagesc(luminance); axis image; colormap(gray);
%%
% To access directly the photons in the image, do this:
photons = sceneGet(sceneMacbethD65,'photons');
%%
% This is a small image, but notice that it is row by col by wavelength
size(photons)
%%
% The values are big because they are photons emitted per second per
% wavelength per steradian per meter from the scene
max(photons(:))
%%
% These are the wavelength sample values in nanometers
wave = sceneGet(sceneMacbethD65,'wave');
%%
% Suppose we compute the mean number of photons across the entire image
meanPhotons = mean(photons,1);  meanPhotons = mean(meanPhotons,2);
meanPhotons = squeeze(meanPhotons);
%%
% and we plot the mean radiance
vcNewGraphWin;
plot(wave,meanPhotons);
xlabel('Wavelength (nm)'); ylabel('Radiance (q/sec/nm/sr/m^2'); 
grid on
%% sceneDescriptions
% To see a general description of the scene, the one printed in the upper
% right of the window, we use this
sceneDescription(sceneMacbethD65)
%%
% Various quantites are stored, such as the horizontal field
% of view in degrees
fprintf('FOV: %f\n',sceneGet(sceneMacbethD65,'fov'))
%%
% If you would like to change the field of view, you can type
sceneMacbethD65 = sceneSet(sceneMacbethD65,'fov',20);
fprintf('FOV: %f\n',sceneGet(sceneMacbethD65,'fov'))
%% Different types of scenes
% There are many types of scenes.  Here is a simple one that is useful for
% demosaicing.  For a list, type help sceneCreate
sceneTest = sceneCreate('freqorientpattern');
vcAddAndSelectObject('scene',sceneTest);
sceneWindow;
%%
% With this one, we might try some simple plots, such as a plot of the
% luminance across the bottom row
sz = sceneGet(sceneTest,'size');
plotScene(sceneTest,'luminance hline',sz);
%%
% We can do this ourselves by getting the luminance of this bottom row as
% follows
luminance = sceneGet(sceneTest,'luminance');
data = luminance(sz(1),:);
%
support = sceneSpatialSupport(sceneTest,'mm');
vcNewGraphWin;
plot(support.x,data,'-');
xlabel('mm'); ylabel('cd/m2'); grid on
%%
rows = round(sceneGet(sceneTest,'rows')/2);
plotScene(sceneTest,'radiance hline',[1,rows]);

%% End of Script









