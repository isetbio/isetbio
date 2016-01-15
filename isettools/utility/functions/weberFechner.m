function y = weberFechner(coef, x)
% y = weberFechnerTVI(coef, x)
%
% Compute the function
% 
%   y = log10(1 ./ (1 + x / coef(1)));
%
% 1/8/16  dhb  Added comment.

% Do it.
y = log10(1 ./ (1 + x ./ coef(1)));
