function r = ieExprnd(mu,varargin)
%Random arrays from exponential distribution.
%
%   R = ieExprnd(MU) returns an array of random numbers chosen from the
%   exponential distribution with location parameter MU.  The size of R is
%   the size of MU.
%
%   R = EXPRND(MU,[M,N,...]) returns an M-by-N-by-... array.
%
%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.
%
% Example:
%   r = ieExprnd(10,1,5000);
%   vcNewGraphWin; hist(r,50);
%
% Modified from Fred Rieke's version of the Mathworks code.
% http://rieke-server.physiol.washington.edu/People/Fred/Classes/545/matlab/StochasticProcessesTutorial/exprnd.m
% 

if nargin < 1
    error('stats:exprnd:TooFewInputs','Requires at least one input argument.');
end

% [err, sizeOut] = statsizechk(1,mu,varargin{:});
% if err > 0
%     error('stats:exprnd:InputSizeMismatch','Size information is inconsistent.');
% end

% Return NaN for elements corresponding to illegal parameter values.
if mu < 0, error('Exponential parameter less than zero.'); end
% mu(mu < 0) = NaN;

sz = zeros(1,length(varargin));
for ii=1:length(varargin);
    sz(ii) = varargin{ii};
end

r = zeros(sz);
% Generate uniform random values, and apply the exponential inverse CDF.
r = -mu .* log(rand(sz)); % == expinv(u, mu)

end
