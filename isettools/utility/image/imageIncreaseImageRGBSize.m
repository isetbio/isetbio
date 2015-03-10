function t = imageIncreaseImageRGBSize(im,s)
% Increase the size of an rgb-style image (r,c,w) by pixel replication.
%
%    t = imageIncreaseImageRGBSize(im,s)
%
% The parameter s is the scale factor.  If the input image is [1,1,w], the
% output image is [s, s, w].   
%  
% Example:
%    [img,p] = imageHarmonic;
%    t = imageIncreaseImageRGBSize(img,3);
%    imagesc(t); colormap(gray)
%
% Copyright ImagEval Consultants, LLC, 2003.

assert(ndims(im) == 3 || ismatrix(im), 'Unexpected im dimension');
if isscalar(s), s = round([s, s]); end

[r,c,~] = size(im);

t = zeros(r*s(1), c*s(2), size(im, 3));
for ii=1:size(im,3)
    t(:,:,ii) = kron(im(:,:,ii),ones(s(1),s(2))); 
end

end