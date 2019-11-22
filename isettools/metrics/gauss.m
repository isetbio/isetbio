function gauss = gauss(hwhm, support)
% Gaussian vector with halfwidth half max and support specified
%
% Syntax:
%   g = gauss(hwhm, support)
%
% Description:
%    The hwhm must be greater than one.
%
%    The hwhm specifies the support of the gaussian between the points
%    where it obtains half of its maximum value. The support indicates the
%    gaussians support in pixels.
%
%    The univariate Gaussian is:
%
%       g1 = (1/s sqrt(2pi)) exp(-(x/2s)^2)
%
%    The univariate hwhm, h, is the value where the Gaussian is at half
%    of its maximum value.   
%
%    The support indicates the gaussians spatial support in pixels. The
%    hwhm must be greater than one. 
%
%    The relationship between the standard deviation, s, & the half max is:
%
%        s  = h / (2*sqrt(ln(2))),  for 1D Gaussian and
%
% Inputs:
%    hwhm    - Numeric. The half width half max point of the gaussian. Note
%              that this value must be > 1.
%    support - Numeric. The spatial support of the gaussian in pixels.
%
% Outputs:
%    gauss   - Matrix. A vector matrix containing the Gaussian vector.
%
% Optional key/value pairs:
%    None.
%

if (nargin < 2),  error('Two input arguments required'); end
if hwhm <= 1, error('Half width half max must be greater than 1.'); end

x = (1:support) - round(support / 2);
s = ieHwhm2SD(hwhm, 1);
gauss = exp(-(x / (2 .* s)) .^ 2);
gauss = gauss / sum(sum(gauss));

% return;
end