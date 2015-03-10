% s_opticsRTPSFView
%
% Create a shift-variant point spread function (svPSF).  Examine the PSFs
% stored in the precomputed svPSF.
%
% (c) Imageval Consulting, LLC, 2012

%%
s_initISET
wbStatus = ieSessionGet('waitbar');
ieSessionSet('waitbar','on');

%% This is a pretty fast way to build a shift variant PSF to explore

oi = oiCreate;
rtOptics = []; spreadLimits = [1 5]; xyRatio = 1.6;
rtOptics = rtSynthetic(oi,rtOptics,spreadLimits,xyRatio);
oi = oiSet(oi,'optics',rtOptics); 

% Get the shift variant PSFs that were computed and stored
angStepSize = 10;
svPSF = rtPrecomputePSF(oi,angStepSize);

%% Alternatively, you can run an oiCompute in the ray trace mode
% The svPSF structure will be stored in the oi. You can get the shift
% variant PSFs that were computed and stored this way.

scene = sceneCreate('point array',384);
scene = sceneSet(scene,'h fov',4);
scene = sceneInterpolateW(scene,550:100:650);
vcAddAndSelectObject(scene); sceneWindow;

oi = oiCompute(oi,scene);
svPSF = oiGet(oi,'psf struct');  
%
%% These are PSFs at the various field heights
vcNewGraphWin; 
for ii=1:size(svPSF.psf,2), imagesc(svPSF.psf{1,ii,1}), axis image; pause(0.3); end

%% These are PSFs at the various sample angles
vcNewGraphWin; 
for ii=1:size(svPSF.psf,1), imagesc(svPSF.psf{ii,end,1}); axis image; pause(0.3); end

%%  Here is a mesh of a point spread at two field heights

% At the center
rtPlot(oi,'psf',550,0);

% At the largest field height
imgHeight = oiGet(oi,'psf image heights','mm');
rtPlot(oi,'psf',550,max(imgHeight(end)));

%%
ieSessionSet('waitbar',wbStatus);

%%