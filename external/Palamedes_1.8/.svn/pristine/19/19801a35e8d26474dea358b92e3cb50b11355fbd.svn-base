%
%PAL_logQuick   Evaluation of the 'Quick' Psychometric Function when 
%   intensity is expressed in log units. (Quick, Kybernetic, 16, pp. 65-67,
%   1974). PAL_logQuick is identical to PAL_Gumbel except that base 2 is 
%   used instead of the natural base.
%
%   syntax: y = PAL_logQuick(params, x)
%
%   y = PAL_logQuick(params, x), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = PAL_logQuick(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y.
%
%   dydx = PAL_logQuick(params, x, 'Derivative') returns the
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
%   y = PAL_logQuick([1 2 0 0], 1)
%
%   y = 0.5
%
%Introduced: Palamedes version 1.6.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function y = PAL_logQuick(params,x,varargin)

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

if ~isempty(varargin)
    if strncmpi(varargin{1}, 'Inverse',3)
        c = (x-gamma)./(1 - gamma - lambda) - 1;
        c = -1.*log(-1.*c)./log(2);
        c = log10(c);
        y = alpha + c./beta;    
    end
    if strncmpi(varargin{1}, 'Derivative',3)
        y = (1 - gamma - lambda).*2.^(-1.*10.^(beta.*(x-alpha))).*log(10).*10.^(beta.*(x-alpha)).*beta;
    end
else
    y = gamma+(1 - gamma - lambda).*(1-2.^(-1.*10.^(beta.*(x-alpha))));
end