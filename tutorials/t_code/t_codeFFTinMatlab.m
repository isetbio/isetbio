%% t_codeFFTinMatlab.m
%
% Explains the Matlab conventions for transforming from space to the frequency domain.

%% Initialize
ieInit;

%% First an extremely small example
nSamples = 6;

%% Inverse transform
%
% In the transform domain, t(1,1) represents the DC term.  You can prove
% this by calculating the inverse FFT for all zeros except t(1,1)
t = zeros(nSamples,nSamples);
t(1,1) = 1;
ft = ifft2(t);
isreal(ft)  % The entries are all 1/(6*6)

%% FFT
%
% In the space domain, the s(1,1) position represents the center of the
% image.  You can prove this by calculation
s = zeros(nSamples,nSamples);
s(1,1) = 1;
fft2(s)     % This produces the output for an impulse at the center
isreal(s)

% The implications of these representations for using fft2 
%
% See Matlab documentation on fft2, ifft2, fftshift and ifftshift

% [DHB Note]: I don't understand the initial comments about centering here.
% I think the center is at floor(N/2)+1 for both odd and even dimension.
% This makes no difference for even, but does change things for odd. Isn't
% the center at [3,3] for N = 5?
%
% Image centering
% The center of an image of size (N,N) is
%  when N is odd,         ceil(N/2) + 1, 
%      so if N = 5, center is (4,4)
%  when N is even is also ceil(N/2) + 1 = N/2 + 1, 
%      so if N = 6, center is (4,4)
%
% When we pad an image or a filter, we want to do so in a way that the
% value at the center remains the center.   
%
% This leaves the center of the image at the same location as it was. So,
% if we have an odd size image, we pad it on the right and bottom. 
%
% Suppose we have an even size image, say N=6, and we pad it to N=7. The
% new center will be at (5,5). To preserve the old center data, which was
% at (4,4) to be at (5,5), we must pad at the left and top first.
%
%  If we pad an odd dimension, first pad the right and bottom
%  If we pad an even dimension, first pad the left and top
%
% These commented out lines might be examples of ways to do things.
%   img = ones(64,64); img = padarray(img,[32 32]);
%   g = zeros(128,128); g(65,65) = 1;

%% PSF/OTF example
%
% Suppose we create a PSF.  In most coding, the natural way to create a
% PSF is as an image.  The center is not in (1,1), but in the center. 
theDim = 128;
g = fspecial('gaussian',theDim,2);
figure(1); subplot(1,3,1); colormap(gray); mesh(g);

% [DHB Note]: I think this should be ifftshift here.  This doesn't matter
% if the support has even dimension, but will matter if it has odd
% dimension.  (fftshift moves the corner to the center, ifftshift moves the
% center to the corner, as I understand it.)  That said, reversing the
% order of the calls does not seem to have any effect.  Note for example
% that in opticsGet, the code to return the center-centered psf from the
% corner centered otf is: val = fftshift(ifft2(otf)); That is, fftshift is
% applied to move the corner to the center.
%
% To calculate the OTF of the point spread function, we should place the
% center of the image in the (1,1) position.  We do this using fftshift.
% We can then take the fft2 of the result to produce the OTF.
gFT = fft2(fftshift(g));
subplot(1,3,2); mesh(abs(gFT)); 

% And now show that you can go back to the original image
gFTAndBack = ifftshift(ifft2(gFT)); 
subplot(1,3,3); mesh(abs(gFTAndBack)); 

%% Image example
% Again, the image center is not in (1,1).  It is in the center.
tmp = load('trees');
cmap = gray(128);
img = cmap(tmp.X);
img = img(1:theDim,1:theDim);
figure(2); subplot(1,4,1); colormap(gray); imagesc(img); axis image

% [DHB Note]: Again, I think this should be ifftshift.
%
% Before we transform the image, we want to place its center in the (1,1)
% position.
imgC = fftshift(img);
subplot(1,4,2); imagesc(imgC); axis image

% Then we compute the transform
imgFT = fft2(imgC);

% We are ready to multiply the transformed image and the OTF 
imgFTgFT = imgFT .* gFT;

% We can return the transform to the space domain.  
imgConvG = ifft2(imgFTgFT);

% When we do, the image center is in the (1,1) position.
subplot(1,4,3); colormap(gray); imagesc(imgConvG); axis image

% [DHB Note]:  I think this should be fftshift, not ifftshift.
%
% We want the center in the center.  So we apply the ifftshift.
imgConvGCentered = ifftshift(imgConvG);
subplot(1,4,4); imagesc(imgConvGCentered); axis image

%% End
