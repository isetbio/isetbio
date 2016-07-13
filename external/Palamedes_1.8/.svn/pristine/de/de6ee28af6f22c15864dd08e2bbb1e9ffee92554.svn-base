%
%PAL_Gumbel   Evaluation of Gumbel Psychometric Function (Gumbel is
%   commonly referred to as 'log Weibull' or, confusingly, as 'Weibull').
%
%   syntax: y = PAL_Gumbel(params, x)
%
%   y = PAL_Gumbel(params, x), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = PAL_Gumbel(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y.
%
%   dydx = PAL_Gumbel(params, x, 'Derivative') returns the
%   derivative (slope of tangent line) of the Psychometric Function
%   evaluated at x.
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.
%
%   A lot of confusion exists regarding the Weibull function family, 
%   please visit: www.palamedestoolbox.org/weibullandfriends.html 
%   if you are not sure whether this is the function you are looking for.
%
%   This example returns the function value at threshold when gamma 
%   ('guess-rate') and lambda ('lapse-rate') both equal 0:
%       
%   y = PAL_Gumbel([1 2 0 0], 1)
%
%   y = 0.6321
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.2.0, 1.6.3 (see History.m)

function y = PAL_Gumbel(params,x,varargin)

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

if ~isempty(varargin)
    if strncmpi(varargin{1}, 'Inverse',3)
        c = (x-gamma)./(1 - gamma - lambda) - 1;
        c = -1.*log(-1.*c);
        c = log10(c);
        y = alpha + c./beta;    
    end
    if strncmpi(varargin{1}, 'Derivative',3)
        y = (1 - gamma - lambda).*exp(-1.*10.^(beta.*(x-alpha))).*log(10).*10.^(beta.*(x-alpha)).*beta;
    end
else
    y = gamma+(1 - gamma - lambda).*(1-exp(-1.*10.^(beta.*(x-alpha))));
end