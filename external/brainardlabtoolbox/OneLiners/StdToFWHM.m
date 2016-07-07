function fullWidthHalfMax = StdToFWHM(standardDeviation)
% fullWidthHalfMax = StdToFWHM(standardDeviation)
% 
% Convert the standard deviation of a Gaussian to its full width at
% half max.
%
% See http://mathworld.wolfram.com/GaussianFunction.html for
% formula.
%
% 4/12/12  dhb  Made this its own function

fullWidthHalfMax = standardDeviation * (2*sqrt(2*log(2)));
