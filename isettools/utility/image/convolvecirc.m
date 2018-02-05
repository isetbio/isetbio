function Y = convolvecirc(X, h)
% Performs 2D circular convolution 
%
% Syntax:
%	Y = convolvecirc(X, h) 
%
% Description:
%    The matrix h (kernel) is convolved with the matrix X. The result has
%    the same size as X. There is probably a Matlab circular convolution by
%    now in the image processing toolbox.
%
%    It is assumed that both the row dimension and column dimension  of h
%    do not exceed those of X.  The result is the same as zero-padding h
%    out to the size of X, and then computing the convolution X * h.
%
% Inputs:
%	 X - The matrix to convolve
%    h - The kernel matrix
%
% Outputs:
%    Y - The resulting matrix, of a size with X
%
% Optional key/value pairs:
%    None.
%

% History:
%	 xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/15/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

[m, n] = size(X);
Y = conv2(X, h);

[r, s] = size(Y);
Y(1:(r - m), :) = Y(1:(r - m), :) + Y((m + 1):r, :);
Y(:, 1:(s - n)) = Y(:, 1:(s - n)) + Y(:, (n + 1):s);

% Clip the extent of the result
Y = Y(1:m, 1:n);

end