function y = ieNormpdf(x,mu,sigma)
%ieNormpdf Normal probability density function (pdf).
%   Y = ieNormpdf(X,MU,SIGMA) Returns the normal pdf with mean, MU, 
%   and standard deviation, SIGMA, at the values in X. 
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.26.
%
%   Examples:
%
%    ieNormpdf(0)
%    x = linspace(-3,3,100); mu = zeros(size(x)); sigma = ones(size(x));
%    n = ieNormpdf(x,mu,sigma); vcNewGraphWin; plot(x,n); grid on
%
%    sigma = 0.2*sigma; mu = mu + 2; x = x + 2;
%    n = ieNormpdf(x,mu,sigma); vcNewGraphWin; plot(x,n); grid on
%
% JRG/BW ISETBIO Team
%

% By default a standard normal
if nargin < 3, sigma = 1; end
if nargin < 2; mu = 0; end
if nargin < 1, 
    error('Requires at least one input argument.');
end

% Should do a size check that arguments are all the same
% if errorcode > 0
%     error('Requires non-scalar arguments to match in size.');
% end

if ~isempty(find(sigma <= 0, 1)), error('Sigma must be > 0'); end

xn = (x - mu) ./ sigma;
y = exp(-0.5 * xn .^2) ./ (sqrt(2*pi) .* sigma);

end