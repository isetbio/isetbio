function g = sumGauss(params, dimension)
% Calculated a weighted sum of three Gaussians, used in SCIELAB
%
% Syntax:
%   g = sumGauss(params, dimension)
%
% Description:
%    Calculate a weighted sum of three gaussians, used in SCIELAB.
%
% Inputs:
%    params    - Array. An array of the required numerical values to
%                calculate the gaussian sum: support, halfwidth1, weight1,
%                halfwidth2, weight2, halfwidth3, and weight3.
%    dimension - Numeric. A number to indicate whether required sum of
%                gaussians is 1-D or 2-D. (Wher 1 = 1D, 2 = 2D) Default 1.
%
% Outputs:
%    g         - Matrix. A matrix representing the weighted gaussian sum.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   scPrepareFilters, gauss, gauss2
%

% History:
%    XX/XX/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/20/19  JNM  Documentation pass

%{
    support = 64;
    halfwidth1 = 5;
    weight1 = 0.2;
    halfwidth2 = 10;
    weight2 = 0.2;
    halfwidth3 = 20;
    weight3 = 0.5;
    params = [support, halfwidth1, weight1, halfwidth2, weight2, ...
        halfwidth3, weight3];
    dimension = 2;
    g = sumGauss(params, dimension);
    X = (1:support) - support / 2;
    vcNewGraphWin;
    mesh(X, X, g);
    grid on;
%}

if notDefined('params'), error('params required'); end
if notDefined('dimension'), dimension = 1; end

width = ceil(params(1));
nGauss = (length(params) - 1) / 2;

if dimension == 2, g = zeros(width, width); else, g = zeros(1, width); end

% Add up the Gaussians.  These calls could be made using fspecial and the
% ieHwhm2SD function as well.  But we leave these here for now because in
% the future we might use oriented Gaussians.
for ii = 1:nGauss
  halfwidth = params(2 * ii);
  weight = params(2 * ii + 1);
  if dimension == 2
    g0 = gauss2(halfwidth, width, halfwidth, width);
  else
    g0 = gauss(halfwidth, width);
  end
  g = g + weight * g0;
end

% Make sure the summed Gaussians equals precisely to 1.
g = g / sum(g(:));

return;