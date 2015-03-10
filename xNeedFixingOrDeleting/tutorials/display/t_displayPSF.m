%% t_displayPSF
%
% HJ is updating the display structure calls to account for display point
% spread functions.  The purpose of the updated structures is to allow us
% to have more spatially accurate descriptions of the display spectral
% radiance.
%
% This script is a includes material about the PSF calculations in the new
% display structure.
%
% NOTES:
%   1) Broken becuase ISETBIO's displayGet does not (yet) have a psf method.
%
% (BW) May 2014

%% Init
s_initISET

%% Make a harmonic image for a display with a psf
I = 0.5*(sin(2*pi*(1:32)/32)+1); I = repmat(I,32,1);
[outImage, d] = displayCompute('LCD-Apple', I);
vcNewGraphWin; imagescRGB(outImage);

%%  Show the display psfs
psf = displayGet(d,'psf');
vcNewGraphWin([],'tall');
for ii=1:3, subplot(3,1,ii), mesh(psf(:,:,ii)); end

%%  Now, replace the psfs with a different shape
for ii=1:3
    psf(:,:,ii) = psf(:,:,ii)';
end
d2 = displaySet(d,'psf',psf);

psf = displayGet(d2,'psf');
vcNewGraphWin([],'tall');
for ii=1:3, subplot(3,1,ii), mesh(psf(:,:,ii)); end

%%  Display the same scene, but with a different psf for the display
outImage = displayCompute(d2, I);
vcNewGraphWin; imagescRGB(outImage);

%% Now do a similar calculation using some of the ISET tools
imSize = [32 32];
scene = sceneCreate('slanted bar',imSize);
I = sceneGet(scene,'rgb');
vcNewGraphWin; imagescRGB(displayCompute(d, I));
vcNewGraphWin; imagescRGB(displayCompute(d2, I));

%% End