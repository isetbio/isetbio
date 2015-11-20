function [img, parms] = imageHarmonicColor(parms)
%Create a windowed spatial harmonic image 
%
%   [img,parms]  = imgHarmonic(parms)
%
% Creates a windowed, oriented, spatial harmonic
%
%      contrast*window.*cos(2*pi*f*([cos(ang)*X + sin(ang)*Y] + ph)) + 1
%
% The parameters are  in the structure parms. See the fields in the
% example, below. 
%
% When some of the fields are vectors (freq, contrast, ang, phase) then the
% return produces the sum of these harmonics.  The sum always has a mean of
% 1.  
%
% The Gabor Flag is used to set the window values (a Gaussian).
% When the flag is non-zero, the value specifies the standard deviation
% of the Gaussian as a fraction of the image size.  For example, if the
% image size is 128 and GaborFlag = 0.25, the standard deviation is 32.
%
% Default parameters are applied if parms is not sent in.  You can see the
% defaults by requesting them on return as below.
%
% Example:
%   [img,p] = imageHarmonic;
%   figure; imagesc(img), colormap(gray); axis image
%
%   parms.row = 32; parms.col = 32; parms.contrast = 1; 
%   parms.ph = pi/2; parms.freq = 2; parms.ang = pi/6;
%   parms.GaborFlag = 0.2;
%   [img,p] = imageHarmonic(parms);
%   vcNewGraphWin; imagesc(img), colormap(gray); axis image
%
% Now, for a sum of two harmonics
%   parms.freq(2) = 3; parms.ang(2) = parms.ang(1);
%   parms.contrast(2) = 1; parms.ph(2) = pi/2;
%   [img,p] = imageHarmonic(parms);
%   vcNewGraphWin; imagesc(img), colormap(gray); axis image
%   plot(img(16,:))
% 
% Copyright ImagEval Consultants, LLC, 2003.

% parms = paramsGaborColorOpponent()

if ~exist('parms','var'), parms = []; end

try ang = parms.ang; catch, ang = 0; parms.ang = ang; end
try contrast = parms.contrast; 
catch, contrast = 1; parms.contrast = contrast; end
try freq = parms.freq; catch, freq = 1; parms.freq = freq; end
try ph = parms.ph; catch, ph = pi/2; parms.ph = ph; end
try row = parms.row; catch, row = 64; parms.row = row; end
try col = parms.col; catch, col = 64; parms.col = col; end
try color = parms.color; catch, color = 1; parms.color = 1; end

% The Gabor Flag is a non-zero value that specifies the standard deviation
% of the Gaussian as a fraction of the image size.  For example, if the
% image size is 128 and GaborFlag = 0.5, the standard deviation is 64.
if isfield(parms,'GaborFlag')
    GaborFlag = parms.GaborFlag; 
else
    GaborFlag = 0; 
    parms.GaborFlag = GaborFlag; 
end

% Calculate the harmonic
[X,Y] = meshgrid((0:(col-1))/col,(0:(row-1))/row);

% Calculate the gabor window
if GaborFlag
    hsize = size(X);
    sigma = GaborFlag*min(row,col);
    g = fspecial('gauss',hsize,sigma);
    g = g/max(g(:));
else 
    g = ones(size(X));
end

% Select the color opponent stimuli

switch parms.color
    case {'siso',1}
        colorOpponent = [0.05 0.10 0.99];
    case{'lmred',2}
        colorOpponent = [0.86 -0.50 -0.05];
    case{'lmsorange',3}
        colorOpponent = [0.12 -0.19 -0.97];
    case{'lmsmagenta',4}
        colorOpponent = [0.12 -0.19 0.97];
end

% harmonics are (1 + sum(cos(2*pi*f*x + ph))
% with the additional issue that X and Y can be at some angle.
img = zeros([size(X) 3]);
for rgbIndex = 1:3
for ii=1:length(freq)
    img(:,:,rgbIndex) = img(:,:,rgbIndex) + ...
        colorOpponent(rgbIndex)*contrast(ii)*g.*cos( 2*pi*freq(ii)*(cos(ang(ii))*X + sin(ang(ii))*Y) ...
        + ph(ii)) + 1;
end
end
img = 0.5 * img / length(freq);


% scfac = sqrt(sum(reshape(img.^2,size(img,1)*size(img,2),3),2));
% scfacrs = reshape(scfac,size(img,1),size(img,2));
% scfacrsrm = repmat(scfacrs,1,1,3);
% img = (sqrt(3)/2)*img./scfacrsrm;

if min(img(:) < 0)
    warning('Harmonics have negative sum, not realizable');
end

% figure; imagesc(img); colormap(gray); axis image

end