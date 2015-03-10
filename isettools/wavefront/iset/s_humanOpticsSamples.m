% Script:  s_humanOpticsSamples
%
% Use the WVF methods to create a human optics in ISET
%
%

%%
s_initISET
scene = sceneCreate('gridlines');
scene = sceneSet(scene,'fov',1);

%% Create wavefront samp_mean human
pupilMM = 4.5; zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wave = (400:10:700)';
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
% wvfP = wvfComputePSF(wvfP);

% Computes the PSF
[d, wvfP] = wvf2PSF(wvfP);

% Look at the OTFs at different wavelengths
% wvfPlot(wvfP,'otf','mm',1,100);
% wvfPlot(wvfP,'otf','mm',2,400);
% wvfPlot(wvfP,'otf','mm',3,300);

oiW = oiCreate('human');
wOptics = siSynthetic('custom',oiW,d);
oiW = oiSet(oiW,'optics',wOptics);
% vcNewGraphWin; plotOTF(oiW,'psf',450);

% vcNewGraphWin; wvfPlot(wvfP,'2d psf space','um',2,50);
% set(gca,'xtick',[-40:10:40],'ytick',[-40:20:40]);
% set(gca,'xlim',[-50 50],'ylim',[-50 50]);

% wvfPlot(wvfP,'otf','mm',1,[]);

%% The version from Marimont
oiM = oiCreate('human');

% vcNewGraphWin; plotOTF(oiM,'otf',550);
% wvfPlot(wvfP,'otf','mm',2,400);

% vcNewGraphWin; plotOI(oiM,'psf',[],550);
% set(gca,'xtick',[-40:10:40],'ytick',[-40:20:40]);
% set(gca,'xlim',[-50 50],'ylim',[-50 50]);

oiM = oiCompute(scene,oiM); oiM = oiSet(oiM,'name','Marimont');
vcAddAndSelectObject(oiM); oiWindow;

oiW = oiCompute(scene,oiW); oiW = oiSet(oiW,'name','Wavefront');
vcAddAndSelectObject(oiW); oiWindow;


