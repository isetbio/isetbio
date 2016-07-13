%
%PAL_HyperbolicSecant   Evaluation of Hyperbolic Secant Psychometric
%   Function
%
%   syntax: y = PAL_HyperbolicSecant(params, x)
%
%   y = PAL_HyperbolicSecant(params, x), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = PAL_HyperbolicSecant(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y.
%
%   dydx = PAL_HyperbolicSecant(params, x, 'Derivative') returns the
%   derivative (slope of tangent line) of the Psychometric Function
%   evaluated at x.
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.
%
%   This example returns the function value at threshold when gamma 
%   ('guess-rate') and lambda ('lapse-rate') both equal 0:
%       
%   y = PAL_HyperbolicSecant([1 2 0 0], 1)
%
%   y = 0.5000
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.1, 1.2.0, 1.6.3 (see History.m)

function y = PAL_HyperbolicSecant(params,x,varargin)

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

if ~isempty(varargin)
    if strncmpi(varargin{1}, 'Inverse',3)
        c = (x - gamma)./(1 - gamma - lambda);
        y = alpha + 2.*log(tan(pi.*c./2))./(pi.*beta);
    end
    if strncmpi(varargin{1}, 'Derivative',3)
        y = (1 - gamma - lambda).* (1+(exp((pi/2).*beta.*(x-alpha))).^2).^-1.*exp((pi/2).*beta.*(x-alpha)).*beta;
    end
else
    y = gamma + (1 - gamma - lambda).*((2/pi).*atan(exp((pi/2).*beta.*(x-alpha))));
end