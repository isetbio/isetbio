% Walk through the calculations in oiCompute.
%
% Description:
%    Walk through the calculations in oiCompute.
%
%    This illustrates how scene radiance is converted through a lens to an
%    optical image (irradiance)
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010.
%    12/21/17  dhb  Clear fSupport and sSupport above line 75, so this
%                   runs. It was broken when I got to it.
%    12/21/17  dhb  Ablate direct calls to fft2/ifft2 in deference to
%                   common routine.
%    09/10/18  jnm  Formatting

%% This is the basic radiance to irradiance code 
% Creates an array of points
scene = sceneCreate('point array');
scene = sceneSet(scene, 'hfov', 0.5);
ieAddObject(scene);
sceneWindow;

% Diffraction limited optics
oi = oiCreate;

% Compute optical image and show it
oi = oiCompute(scene, oi);
ieAddObject(oi);
oiWindow;

%% Make a bigger f-number, compute and show
optics = oiGet(oi, 'optics');
fnSmall = opticsGet(optics, 'f number');
fnBig = 3 * fnSmall;

optics = oiGet(oi, 'optics');
optics = opticsSet(optics, 'fNumber', fnBig);
oi2 = oiSet(oi, 'optics', optics);
oi2 = oiCompute(scene, oi2);
ieAddObject(oi2);
oiWindow;

%% Plot the psf of the optics
vcNewGraphWin;
thisWave = 600;
oiPlot(oi, 'psf', [], thisWave);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20]);
colormap(0.5 * copper + 0.5 * ones(size(copper)));

%% Plot irradiance image
vcNewGraphWin;
gridSpacing = 5;  % um
oiPlot(oi, 'irradiance image with grid', [], gridSpacing);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20])
title(sprintf('F-number = %d', fnSmall))

%% What happens if we change the f/# of the optics and replot?
vcNewGraphWin;
oiPlot(oi2, 'psf', [], thisWave);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20])

colormap(0.5 * copper + 0.5 * ones(size(copper)))
title(sprintf('F-number = %d', fnBig))

%% Plot new irradiance image
vcNewGraphWin;
gridSpacing = 5;
oiPlot(oi2, 'irradiance image with grid', [], gridSpacing);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20])
title(sprintf('F-number = %d', fnBig))

%% Here is the psf plot method, including the OTF and PSF
% This is just copied from the oiPlot code, really.
units = 'um';  % Specify units!

% Calling conventions should be specified here. The opticsGet() for
% diffraction limited should be adjusted so that this code becomes shorter.
% idx = ieFindWaveIndex(wavelength, thisWave);
clear fSupport sSupport;
nSamp = 100;  % Number of frequency steps from 0 to incoherent cutoff
val = opticsGet(optics, 'dlFSupport', thisWave, units, nSamp);
[fSupport(:, :, 1), fSupport(:, :, 2)] = meshgrid(val{1}, val{2});

% We over sample the frequency to get a smoother PSF image. You can specify
% the factor for oversampling if you like in the calling arguments.
s = 4;
fSupport = fSupport * s;

% Frequency units are cycles/micron The spatial frequency support runs from
% -Nyquist:Nyquist. With this coding, the Nyquist frequency is actually the
% peak frequency value. There are two samples per Nyquist, so the sample
% spacing is 1/(2*peakF)
deltaSpace = 1 / (2 * max(fSupport(:)));

% Diffraction limited OTF/PSF
otf = dlMTF(oi, fSupport, thisWave, units);
[~, ~, psf] = OtfToPsf([], [], fftshift(otf));

samp = (-nSamp:(nSamp-1));
[X, Y] = meshgrid(samp, samp);
sSupport(:, :, 1) = X * deltaSpace;
sSupport(:, :, 2) = Y * deltaSpace;

fNumber = opticsGet(optics, 'fNumber');

% First zero crossing
radius = (2.44 * fNumber * thisWave * 10 ^ -9) / 2 ...
    * ieUnitScaleFactor(units);

nCircleSamples = 200;
[~, ptsXY]=ieShape('circle', 'nSamp', nCircleSamples, 'radius', radius);
adX = ptsXY(:, 1);
adY = ptsXY(:, 2);
adZ = zeros(size(ptsXY(:, 1)));

x = sSupport(:, :, 1);
y = sSupport(:, :, 2);

% We used to choose a size to show that illustrates the PS down to 1%.
% Now we just show the whole thing.
% [Note: JNM - Based on the plot support provided below, isn't that .1%?]
%{
    sz = selectPlotSupport(psf, 0.001);
    x =  getMiddleMatrix(x, sz);
    y =  getMiddleMatrix(y, sz);
    psf = getMiddleMatrix(psf, sz);
%}
clf
mesh(x, y, psf);

% For the diffraction limited case, we draw the Airy disk
hold on;
plot3(adX, adY, adZ, 'k.');
hold off;
colormap(0.5 * copper + 0.5 * ones(size(copper)))
