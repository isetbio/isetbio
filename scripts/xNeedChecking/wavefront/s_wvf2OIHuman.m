%% s_wvf2OIHuman
%
% Create a wavefront objects from the sample mean human wvf data
%
% We tested the relationship between the wvf psf value and the converted
% psf in ISET using s_wvf2OI.m.  In that case, we used the diffraction
% limited instance. That script gives us some modest confidence in the idea
% that the conversion from wvf to ISET is about right.
%
% Here, we convert the human psf, based on the Thibos data, into ISET. We
% explore the human model from Thibos et al. and Marimont and Wandell.
% This script suggests that there are modest differences between the mean
% Thibos measurement and the Marimont and Wandell calculation.  But, they
% are actually pretty close, to my amazement (BW).
%
% (c) Wavefront Toolbox Team, 2012


%% Initialize
s_initISET
maxUM = 10;
scene = sceneCreate; vcAddObject(scene);

%% Convert WVF human data to ISET
% This cell shows that the conversion for  WVF to ISET work well for one
% wavelength.
%
% The data were collected by Thibos and are described in the wvfLoadThibosVirtualEyes
% function and reference therein.
wave = 550; wave = wave(:);
pupilMM = 3; 

% Load human wvf
% 
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfComputePSF(wvfP);

thisWave = wvfGet(wvfP,'calc wave');
[u,p,f] = wvfPlot(wvfP,'2d psf space','um',thisWave,maxUM);
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

oiD = wvf2oi(wvfP,'human');
oiD = oiSet(oiD,'name','Human 3mm');
vcAddObject(oiD); oiWindow;
uData = oiPlot(oiD,'psf','um',thisWave,maxUM);
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

% Not right ... something needs to be fixed here (BW).
% Interpolate the wvf data onto the spatial grid of the ISET data.  Same as
% above, but written more compactly.
test = interp2(u.x,u.y,u.z,uData.x,uData.y);
test = test/sum(test(:));

vcNewGraphWin([],'tall');
subplot(3,1,1), mesh(uData.x,uData.y,uData.psf)
subplot(3,1,2), mesh(uData.x,uData.y,test)
subplot(3,1,3), mesh(uData.x,uData.y,test - uData.psf)

vcNewGraphWin;
plot(uData.psf(:),test(:),'.');
xlabel('ISET psf values');
ylabel('WVF psf values');

grid on

% The conversion between minutes and 'um' is not right for the typical
% human focal length.  Look into this.
% wvfPlot(wvfP,'2d psf space','min',thisWave);


%% Now, one more wavelength
% This is slightly imperfect.  But it is fairly close.  Not sure what the
% source of the imperfection is, and we must track it down. It could be
% because of the 0.25 micron issue discussed in s_wvf2OI.m (BW). 
wave = 500; wave = wave(:);
pupilMM = 3; 

% Load human wvf
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfComputePSF(wvfP);

thisWave = wvfGet(wvfP,'calc wave');
[u,p,f] = wvfPlot(wvfP,'2d psf space','um',thisWave);
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

oiD = wvf2oi(wvfP,'human');
oiD = oiSet(oiD,'name','Human 3mm');
vcAddAndSelectObject(oiD); oiWindow;
uData = oiPlot(oiD,'psf','um',thisWave);
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

% Interpolate the wvf data onto the spatial grid of the ISET data.  Same as
% above, but written more compactly.
test = interp2(u.x,u.y,u.z,uData.x,uData.y);
test = test/sum(test(:));

vcNewGraphWin([],'tall');
subplot(3,1,1), mesh(uData.x,uData.y,uData.psf)
subplot(3,1,2), mesh(uData.x,uData.y,test)
subplot(3,1,3), mesh(uData.x,uData.y,test - uData.psf)
subplot(3,1,1)
title(sprintf('Infocus %.0f Measured %.0f',wvfGet(wvfP,'measured wl'),wvfGet(wvfP,'calc wavelengths')))

vcNewGraphWin;
plot(uData.psf(:),test(:),'.');
grid on

%% Now compute a few wavelengths
% This looks very odd, especially there is a sharp change around 500 to 480
% nm that seems very unlikely.
wave = 400:10:700; wave = wave(:);
pupilMM = 3; 
thisWave = 550;

% Load human wvf
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfComputePSF(wvfP);

% [u,p,f] = wvfPlot(wvfP,'2d psf space','um',thisWave);
% set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

oiD = wvf2oi(wvfP,'human');
oiD = oiSet(oiD,'name','Human WVF 3mm');

oiMW = oiCreate('human');
oiMW = oiSet(oiMW,'name','Human MW 3mm');

%% Now compare the slanted bar response in the OI 
% These are reasonably close for calculations separated by so many years.

scene = sceneCreate('slanted bar');
scene = sceneSet(scene,'h fov',1);
oiD = oiCompute(oiD,scene);
oiD = oiSet(oiD,'name','Thibos');
vcAddAndSelectObject(oiD); oiWindow;
uT = oiPlot(oiD,'irradiance hline',[1,200]);  % (x,y) format
title('Thibos')

oiMW = oiCompute(oiMW,scene);
oiMW = oiSet(oiMW,'name','MW');
vcAddAndSelectObject(oiMW); oiWindow;
uM = oiPlot(oiMW,'irradiance hline',[1,200]);  % (x,y) format
title('Marimont and Wandell')

%% The local error is pretty substantial
% This might matter for calculation properties of the cone absorptions. We
% will see.
vcNewGraphWin;
mesh(uT.pos,uT.wave,(uT.data - uM.data) ./ uM.data);

%% Show the photon absorptions

sensor = sensorCreate('human');
sensor = sensorSet(sensor,'exp time',0.050);
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'hfov'),scene,oiD);

sensorD = sensorCompute(sensor,oiD);
sensorD = sensorSet(sensorD,'name','Thibos calc');
vcAddAndSelectObject(sensorD); sensorWindow('scale',1);

sensorMW = sensorCompute(sensor,oiMW);
sensorMW = sensorSet(sensorMW,'name','MW calc');
vcAddAndSelectObject(sensorMW); sensorWindow('scale',1);

%% End

