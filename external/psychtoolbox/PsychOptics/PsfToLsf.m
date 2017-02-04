function lsf = PsfToLsf(psf,varargin)
%LSFTOPSF  Convert a line spread function to a point spread function
%    lsf = LSFTOPSF(psf,varargin)
%
%    This works by convolving a horizontal line with a psf and returning the
%    vertical slice through the center.  The spatial support of the lsf is
%    equal to the spatial support of the vertical slice through the passed
%    psf.
%
%    The returned lsf is normalized to have a peak amplitude of 1.
%
%    The psf must be passed with even-dimensioned spatial support.
%
%    See also LSFTOPSF.

% Passed row and colum dimension must be even.
[m,n] = size(psf);
if (mod(m,2) ~= 0)
    error('Number of rows must be even');
end
if (mod(n,2) ~= 0)
    error('Number of cols must be even');
end

% Create a 2D image of a line across the center row.
aLine2D = zeros(size(psf));
aLine2D(m/2+1,:) = 1;

% Convolve with the psf
aLineConvolved = conv2(aLine2D,psf,'same');

% Extract the center column and normalize
lsf = aLineConvolved(:,n/2+1);
lsf = lsf/max(lsf(:));

end