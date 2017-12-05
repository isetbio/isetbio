function MP = ieN2MegaPixel(N, precision)
% Compute megapixel count from N
%
% Syntax:
%   mp = ieN2MegaPixel(n, prec)
%
% Description:
%    Compute the megapixel count from N. The decimal precision defaults to
%    1, but you can ask for more.
%
% Inputs:
%    N         - The value
%    precision - The precision point
%
% Outputs:
%    MP        - The megapixel count
%
% Example:
%{
	N = 1024 * 1024;
    MP = ieN2MegaPixel(N)
%}
%{
    N = 1920 * 1024;
    MP = ieN2MegaPixel(N, 2)
    MP = ieN2MegaPixel(N, 0)
%}

if notDefined('N'), error("N must be defined"); end
if notDefined('precision'), precision = 1; end

MP = round(N * 1e-6 * (10 ^ precision)) / 10 ^ precision;

end