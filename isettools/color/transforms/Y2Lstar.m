function Lstar = Y2Lstar(Y,Yn)
% Convert Y (luminance) to L* (CIELAB)
%
% Syntax:
%   Lstar = Y2Lstar(Y, Yn)
%
% Description:
%    Convert luminance (1931 CIE Y coordinate) into CIELAB/CIELUV L* The
%    luminance of the white point (Yn) is required.
%
% Inputs:
%    Y     - Luminance
%    Yn    - Luminance of the white point
%
% Outputs:
%    Lstar - The lightness
%
% References:
%    Wyszecki and Stiles, others.
%
% Copyright ImagEval Consultants, LLC, 2003.

% Basic formula
T = Y / Yn;
Lstar = 116 * (T .^ (1/3)) - 16;

% Buf if the ratio is small ...
l = (T < .008856);
Lstar(l) = 903.3 * T(l);

end