function vol = GetUnitSphereVol(dim)
% vol = GetUnitSphereVol(dim)
%
% Compute the volume of a unit spehre in dim dimensions.
% Formula form Duda and Hart, p. 24.
%
% 6/29/04   dhb     Wrote it.

% Check dimension
if (dim <= 0)
    error('Passed dimension is not positive');
end

% Dimension is odd
if (rem(dim,2) == 1)
    halfDim = (dim-1)/2;
    vol = (2^dim)*(pi^halfDim)*factorial(halfDim)/factorial(dim);
    
% Dimension is even
elseif (rem(dim,2) == 0)
    halfDim = dim/2;
    vol = (pi^halfDim)/factorial(halfDim);
else
    error('Passed dimension is not an integer');
end