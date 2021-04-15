%% Show a letter on a display imaged at the level of the cone mosaic
%
%
% TODO:
%   Write a routine that takes the outline of the letter and superimposes
%   it on the color image of the cones.
%
%   Implement this for the hex mosaic.
%
% Wandell, 2019
%
% See also
%   sceneCreate('letter', ...), displayCreate, .... t_cones*, t_cmosaic*
%

%%
ieInit;

%% Create a letter on a display

% family, size, dpi
sceneFOV = 3;
font = fontCreate('A', 'Georgia', 10, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',sceneFOV);
scene = sceneCombine(scene,scene,'direction','horizontal');
scene = sceneCombine(scene,scene,'direction','vertical');

% We should pad the scene so the eye movements do not move the scene beyond
% the array

% Here is the scene
sceneWindow(scene);

%% Push the scene through human optics

oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
oiWindow(oi);

%%
hFig = figure(1);
% set(hFig, 'Position', [10 10 1400 1200]);
 
% sv = NicePlot.getSubPlotPosVectors(...
%        'rowsNum', 3, ...
%        'colsNum', 3, ...
%        'heightMargin',  0.07, ...
%        'widthMargin',    0.05, ...
%        'leftMargin',     0.05, ...
%        'rightMargin',    0.03, ...
%        'bottomMargin',   0.06, ...
%        'topMargin',      0.03);
   
%%  Now image it on the cone mosaic
 
% cones = coneMosaic;
% Generate mosaic centered at target eccentricity
mosaicEcc = [0 0];

% There are very many parameters in many routines.  It is hard to find all
% of them. I would like this to look like this 
%{  
   % Many important routines should do this
   cmP = cMosaic('params');   % Return the settable params in a struct

   % We choose which ones to set.  The params should be all lowercase with
   % no spaces.
   cmP.sizedegs = [1 1]*0.6;  % We set a couple of them
   cmP.eccentricitydegs = mosaicEcc;

   cm = cMosaic(cmP);         % We call cMosaic with the struct
%}

fov = sceneGet(scene,'fov');
cm = cMosaic('sizeDegs', [fov,fov], ...
    'eccentricityDegs', mosaicEcc); 
 
noiseFreeExcitationResponse = cm.compute(oi);   %% Error here

vParams = cm.visualize('params');

% Visualize mosaic response
cm.visualize( ...
    'domain', 'degrees', ...
    'activation', noiseFreeExcitationResponse, ...
    'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));


