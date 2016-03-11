function y = hill(coef, x)
% y = hill(coef, x)
%
% Compute hill function
%
% y = 1.0 ./ (1 + (coef(1) ./ x).^coef(2));
%
% 1/8/16  dhb  Added comment.

% Do it
y = 1.0 ./ (1 + (coef(1) ./ x).^coef(2));
