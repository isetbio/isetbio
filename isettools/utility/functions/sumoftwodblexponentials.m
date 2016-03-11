function y = sumoftwodblexponential(coef, x)
% y = sumoftwodblexponential(coef, x)
%
% Compute sum of two double exponential functions, with no mean offset.
% 
%  y = A1*exp(-C1*abs(x) + A2*exp(-C2*abs(x);
%
% A1: coef(1)
% C1: coef(2)
% A2: coef(3)
% C2: coef(4)
%
% 1/8/16  dhb  Added commments
%         dhb  Changed both name and what is computed to be more standard

% Do it
y = coef(1) .* exp(-coef(2)*abs(x)) + coef(3) .* exp(-coef(4)*abs(x));

