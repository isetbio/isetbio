function standardDeviation = FWHMToStd(fullWidthHalfMax)
% standardDeviation = FWHMToStd(fullWidthHalfMax)
% 
% Convert the full width half max of a Gaussian to its standard
% deviation.
%
% See http://mathworld.wolfram.com/GaussianFunction.html for
% formula.
%
% 4/12/12  dhb  Made this its own function

standardDeviation = fullWidthHalfMax / (2*sqrt(2*log(2)));
