%% Show the basic quadrature calculation
%
% This is for a 1D image.
% Tomorrow extend the calculations to a 2D image.
%

%%
ieInit

%% Testing the basic idea

%  Make a signal and two harmonics in quadrature
x = [0:127]/128;
f = 2;
s = sin(2*pi*f*x);
c = cos(2*pi*f*x);
sig = 0.4*square(2*f*pi*x) + 0.5;

eBase = dot(sig,s)^2 + dot(sig,c)^2;

%
ieFigure; plot(x,sig,'k-',x,s,'r-',x,c,'b-');

%% Take the inner product of the signal with each harmonic. 

% Then compute the energy, also known as the amplitude at that
% frequency.

% Shift the signal and recompute
for ii=1:2:10
    eShift = dot(circshift(sig,ii),s)^2 + dot(circshift(sig,ii),c)^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now do the same, but for a 2D image
img    = repmat(sig,[128,1]);
simg   = repmat(s,[128,1]);
cimg   = repmat(c,[128,1]);
ieFigure; imagesc(img); colormap(gray); axis image

eBase = dot(img(:),simg(:))^2 + dot(img(:),cimg(:))^2;
for ii=1:2:10
    thisIMG = circshift(img,ii,2);
    imagesc(thisIMG); colormap(gray); axis image; pause(0.2);
    eShift = dot(thisIMG(:),simg(:))^2 + dot(thisIMG(:),cimg(:))^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now modify the calculation by applying by a Gaussian envelope

% Big difference with a small envelope, and little difference with a big
% envelope, like the full harmonic above.
spread = 32;
g = fspecial('gaussian',[128 128],spread);

ieFigure; imagesc(g); colormap(gray); axis image
gsimg = g .* simg;
gcimg = g .* cimg;

eBase = dot(img(:),gsimg(:))^2 + dot(img(:),gcimg(:))^2;

% Small shifts and there is a constant response
fprintf('\n----\n');
for ii=1:2:spread
    thisIMG = circshift(img,ii,2);
    eShift = dot(thisIMG(:),gsimg(:))^2 + dot(thisIMG(:),gcimg(:))^2;
    fprintf('Difference (percentage): %.6f (step %d)\n',100*(eBase - eShift)/eBase,ii)
end
