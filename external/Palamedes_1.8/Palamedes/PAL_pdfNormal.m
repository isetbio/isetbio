%
%PAL_pdfNormal  Normal probability density
%
%syntax: y = PAL_pdfNormal(x, Mean, SD)
%
%Returns the probability density of the normal distribution with mean 
%   'Mean' and standard deviation 'SD' evaluated at 'x'.
%
%   'x' may be array of any size, 'Mean' and 'SD' should be scalars or
%       arrays in size equal to 'x'.
%
%Example:
%
%   y = PAL_pdfNormal([-1 0 1], 0, 1) returns:
%
%   y =
%
%       0.2420    0.3989    0.2420
%
%Introduced: Palamedes version 1.0.0 (NP)

function y = PAL_pdfNormal(x, Mean, SD)

y = exp(-1.*(x-Mean).^2./(2*SD.^2))./((2*pi).^.5*SD);