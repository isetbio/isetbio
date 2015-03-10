%%  s_wvf2OI
%
%  Create an ISET optical image based on the point spread functions
%  calculated with the wavefront toolbox.
%
%  We are currently trying to figure out whether we have the units right.
%  It is surprising to see that the aberration from the wvf is so extreme.
%  It is much more extreme than the Marimont and Wandell 3mm pupil
%  analyses.  If!! we have it right here.
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize

s_initISET

s = which('s_wvf2OI');
cd(fileparts(s));
maxUM = 10;

%% Convert a diffraction limited PSF
% Recall that we compare the diffraction limited case in
% v_wvfDiffractionPSF.m and found agreement between WVF, ISET, and PTB.
% Here we simply test whether we can convert a single wavelength
% pointspread in WVF into the same pointspread in ISET using the wvf2oi
% function.
wvfP = wvfCreate;
wvfP = wvfComputePSF(wvfP);
thisWave = wvfGet(wvfP,'calc wave');

% The conversion between minutes and 'um' is not right for the typical
% human focal length.  Look into this.
% wvfPlot(wvfP,'2d psf space','min',thisWave);
[u,p] = wvfPlot(wvfP,'2d psf space','um',thisWave,maxUM);

%% Create a dummy scene and the OI structure

% The dummy scene is needed mainly to establish the viewing distance
scene = sceneCreate; vcAddObject(scene);

% Create the corresponding OI structure, add it to the database
oiD = wvf2oi(wvfP,'shift invariant');
oiD = oiSet(oiD,'name','Diffraction limited');
vcAddAndSelectObject(oiD);

% Plot the psf
uData = plotOI(oiD,'psf','um',thisWave);

%% Now, compare the plot in uData and u.  

% They should be the same.
% Looks good on July 12, 2012 (BW).

% Interpolate the wvf data onto the spatial grid of the ISET data.
test = interp2(u.x,u.y,u.z,uData.x,uData.y,'linear',0);

% The interpolation doesn't preserve the sum.  We need to deal with this
% better, as DHB pointed out.
test = test/sum(test(:));

% Here are the two meshes after interpolation to the ISET grid at 0.25
% microns.  I am surprised (worried) because the uData grid is not 0.25
% microns.  It is off by a little bit (BW).
vcNewGraphWin([],'tall');
subplot(3,1,1), mesh(uData.x,uData.y,uData.psf)
subplot(3,1,2), mesh(uData.x,uData.y,test)
subplot(3,1,3), mesh(uData.x,uData.y,test - uData.psf)

vcNewGraphWin;
plot(uData.psf(:),test(:),'o');
grid on

%% Diffraction case at multiple wavelength
% This section tests whether things work OK for a few wavelengths at once.
wave = (500:50:600); wave = wave(:);
wvfP = wvfCreate('wave',wave,'name',sprintf('Diffraction-limited'));
wvfP = wvfComputePSF(wvfP);

% Edit the choice of wavelength a few times to compare.  They are all OK, I
% think (BW).
thisWave = wave(3);

% The conversion between 'um' and 'min' is not right.  Look into this.
% wvfPlot(wvfP,'2d psf space','min',thisWave);
wvfPlot(wvfP,'2d psf space','um',thisWave);

% Convert to optical image
oiD = wvf2oi(wvfP,'shift invariant');

oiD = oiSet(oiD,'name','Diffraction limited');
% vcAddAndSelectObject(oiD); 
% oiWindow;
uData = plotOI(oiD,'psf','um',thisWave);

wvfPlot(wvfP,'2d psf space','um',thisWave,maxUM);

oiD = wvf2oi(wvfP,'human');
oiD = oiSet(oiD,'name','Human 3mm');
% vcAddAndSelectObject(oiD); oiWindow;

% Interpolate the wvf data onto the spatial grid of the ISET data.  Same as
% above, but written more compactly.
test = interp2(u.x,u.y,u.z,uData.x,uData.y,'linear',0);
test = test/sum(test(:));

vcNewGraphWin([],'tall');
subplot(3,1,1), mesh(uData.x,uData.y,uData.psf)
subplot(3,1,2), mesh(uData.x,uData.y,test)
subplot(3,1,3), mesh(uData.x,uData.y,test - uData.psf)

vcNewGraphWin;
plot(uData.psf(:),test(:),'o');
grid on

%% END
