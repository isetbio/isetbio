% Account for display PSFs.
%
% Description:
%    This script is designed to help with the understanding of how to
%    account for the display PSFs.
%
%    HJ is updating the display structure calls to account for display
%    point spread functions. The purpose of the updated structures is to
%    allow us to have more spatially accurate descriptions of the display
%    spectral radiance.
%
%    This script is a includes material about the PSF calculations in the
%    new display structure.
%
% Notes:
%    * Broken becuase ISETBIO's displayGet does not (yet) have a psf method
%    * TODO: Determine if HJ's update is complete.
%    * TODO: Determine if this is still broken or not.
%    * TODO: Determine if the further documentation mentioned below for the
%      PSFs was ever added. If not, add.
%

% History:
%    05/xx/14  BW   Created May 2014
%    09/19/18  jnm  Formatting

%% Init
ieInit

%% Make a harmonic image for a display with a psf
x = 1:32;  % Spatial samples
A = 0.25;  % Harmonic amplitude
M = 0.5;   % Mean level

I = A * sin(2 * pi * x / max(x)) + M;
I = repmat(I, max(x), 1);
[outImage, d] = displayCompute('LCD-Apple', I);
vcNewGraphWin;
imagescRGB(outImage);

%% Show the display point spread functions (psfs)
psfs = displayGet(d, 'dixel intensity map');

vcNewGraphWin([], 'tall');
for ii = 1:3, subplot(3, 1, ii), mesh(psfs(:, :, ii)); end
xlabel('x position');
ylabel('y position')
zlabel('Relative intensity');

% We need to explain the psfs scaling and resampling. We carefully set up
% values to preserve the units. Further documentation is needed. Note that
% the sum is xSize * ySize. When we change the spatial sampling, these
% functions change amplitude as well so that the units work out right.
% sum(sum(psfs(:, :, 1)))


%% Now, replace the psfs with a different shape
for ii = 1:3, psfs(:, :, ii) = psfs(:, :, ii)'; end
d2 = displaySet(d, 'dixel intensity map', psfs);

psfs = displayGet(d2, 'dixel intensity map');
vcNewGraphWin([], 'tall');
for ii = 1:3, subplot(3, 1, ii), mesh(psfs(:, :, ii)); end
xlabel('x position');
ylabel('y position')
zlabel('Relative intensity');

%% Display the same scene, but with a different psf for the display
outImage = displayCompute(d2, I);
vcNewGraphWin;
imagescRGB(outImage);

%% Now do a similar calculation using some of the ISET tools
imSize = [32 32];
scene = sceneCreate('slanted bar', imSize);
I = sceneGet(scene, 'rgb');
vcNewGraphWin([], 'tall');
subplot(2, 1, 1), imagescRGB(displayCompute(d, I));
title('Display V');
subplot(2, 1, 2), imagescRGB(displayCompute(d2, I));
title('Display H')

%% End