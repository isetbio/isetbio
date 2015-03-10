function im = checkerboard(checkPeriod,nCheckPairs)
% Create a checkerboard image
%
%   im = checkerboard(checkPeriod,nCheckPairs)
%
% Purpose:
%   A black and white checkerboard image, suitable for using as part of a
%   test scene (say for optical geometric distortion) is returned.  The
%   checkPeriod is specified in terms of pixels.  The number of pairs of
%   black and white checks (both vertically and horizontally) is also
%   created. 
%
% Example:
%   im = checkerboard(16,8); 
%   imshow(im); colormap(gray);
%   imwrite(im,'checkerboard.jpg','jpeg');
%
% Copyright ImagEval Consultants, LLC, 2005.


if notDefined('checkPeriod'), checkPeriod = 16; end
if notDefined('nCheckPairs'), nCheckPairs =  8; end

basicPattern = kron([0,1; 1 ,0], ones(checkPeriod));
im = repmat(basicPattern,nCheckPairs,nCheckPairs);

end