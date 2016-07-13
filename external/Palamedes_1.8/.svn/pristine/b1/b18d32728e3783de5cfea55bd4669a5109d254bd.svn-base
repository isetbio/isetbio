%
%PAL_Weibull   Evaluation of Weibull Psychometric Function
%
%   syntax: y = PAL_Weibull(params, x)
%
%   y = PAL_Weibull(params, x), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = PAL_Weibull(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y.
%
%   dydx = PAL_Weibull(params, x, 'Derivative') returns the
%   derivative (slope of tangent line) of the Psychometric Function
%   evaluated at x.
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.
%
%   Note that the Weibull function expects x-values on a linear, not
%   logarithmic, scale. The function evaluates to the guess rate at x = 0
%   (i.e., x = 0 indicates an absence of signal). Negative values for x are
%   meaningless. Use PAL_Gumbel when x values are log-transformed.
%
%   A lot of confusion exists regarding the Weibull function family, 
%   please visit: www.palamedestoolbox.org/weibullandfriends.html 
%   if you are not sure whether this is the function you are looking for.
%
%   This example returns the function value at threshold when gamma 
%   ('guess-rate') and lambda ('lapse-rate') both equal 0:
%       
%   y = PAL_Weibull([1 2 0 0], 1)
%
%   y = 0.6321
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.1, 1.2.0, 1.6.3 (see History.m)

function y = PAL_Weibull(params,x,varargin)

if min(x) < 0
    message = 'The Weibull function is not defined for negative stimulus intensity values. ';
    message = [message 'In case your stimulus intensity values are log transformed values, you '];
    message = [message 'should use the Gumbel (or ''log-Weibull'') function (PAL_Gumbel).'];
    message = [message 'Visit www.palamedestoolbox.org/weibullandfriends.html for more '];    
    message = [message 'information'];    
    error('PALAMEDES:negativeXweibull',message);
end

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

if ~isempty(varargin)
    if strncmpi(varargin{1}, 'Inverse',3)
        c = (x - gamma)./(1 - gamma - lambda);
        y = alpha.*(-log(1 - c)).^(1./beta);    
    end
    if strncmpi(varargin{1}, 'Derivative',3)
        y = (1-gamma-lambda).*exp(-1*(x./alpha).^beta).*(x./alpha).^(beta-1).*beta./alpha;
    end
else
    y = gamma + (1 - gamma - lambda).*(1 - exp(-1*(x./alpha).^beta));
end