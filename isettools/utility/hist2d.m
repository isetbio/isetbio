function Hout = hist2d(D, Xn, Yn, Xrange, Yrange)
% Calculate and return a 2D histogram of D
%
% Syntax:
%   Hout = hist2d(D, [Xn], [Yn], [Xrange], [Yrange]);
%
% Description:
%    The function calculates and returns a 2D Histogram of D
%
% Inputs:
%    D      - Must be an array of complex numbers, or a 2x matrix (2 rows
%             or 2 columns)
%    Xn     - (Optional) The number of points on the X-axis for the
%             histogram. Default is 20. 
%    Yn     - (Optional) The number of points on the Y-axis for the
%             histogram. Default is 20. 
%    Xrange - (Optional) The range of X data, of the format [xLow xHigh].
%             Default falls to the min and max of the provided data.
%    Yrange - (Optional) The range of Y data, of the format [yLow yHigh].
%             Default falls to the min and max of the provided data.
%
% Outputs:
%    Hout   - The histograph data
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    linspace
%

% Examples:
%{
    hist2d([randn(1, 10000);
    randn(1, 10000)])
%}

% first supply optional arguments
if nargin < 3, Yn = 20; end
if nargin < 2, Xn = 20; end
if ~isreal(D), D = [real(D(:)) imag(D(:))]; end
if (size(D, 1) < size(D, 2) && size(D, 1) > 1), D = D.'; end
if size(D, 2) ~= 2
    error('The input data matrix must have 2 rows or 2 columns');
end
if nargin < 4, Xrange = [min(D(:, 1)), max(D(:, 1))]; end
if nargin < 5, Yrange = [min(D(:, 2)), max(D(:, 2))]; end
%
Xlo = Xrange(1);
Xhi = Xrange(2);
Ylo = Yrange(1);
Yhi = Yrange(2);
X = linspace(Xlo, Xhi, Xn)';
Y = linspace(Ylo, Yhi, Yn)';

Dx = D(:, 1);
Dy = D(:, 2);
n = length(D);

H = zeros(Yn, Xn);

for i = 1:n
    x = dsearchn(X, Dx(i));
    y = dsearchn(Y, Dy(i));
    H(y, x) = H(y, x) + 1;
end

figure, surf(X, Y, H);
% Xmid = 0.5 * (X(1:end - 1) + X(2:end));
% Ymid = 0.5 * (Y(1:end - 1) + Y(2:end));
% figure, pcolor(Xmid, Ymid, H); 

colorbar;
shading flat;
axis square tight; 
if nargout > 0, Hout = H; end

end