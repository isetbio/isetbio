%
%PAL_spreadPF      Derives the spread (or support) of Psychometric Function
%   (PF)
%
%Syntax: spread = PAL_spreadPF(params, delta, PFname)
%
%   The spread of a PF is a measure of the PFs width. Specifically, it is 
%   the range of stimulus intensities (or, more generally, the IV) across 
%   which the PF evaluates from gamma (PF's 'guess-rate') + 'delta' to 
%   1 - lambda (PF's 'lapse-rate') - delta. 'delta' must have value greater 
%   than 0 and less than (1 - gamma - lambda)/2.
%
%   'PFname's supported are: 'Logistic', 'CumulativeNormal', 'Weibull',
%   'Quick', 'Gumbel', 'logQuick', and 'HyperbolicSecant'.
%
%Example:
%   spread = PAL_spreadPF([0 1 0 0], (1 - .6826)/2, 'CumulativeNormal')
%
%   returns:
%
%   spread = 1.9996
%
%   (confirming the well-known fact that approximately 68.26% of scores in 
%       the normal distribution falls between Z = -1 and Z = 1).
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.4.0, 1.6.1, 1.6.3 (see History.m)

function spread = PAL_spreadPF(params, delta, PFname)

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

switch lower(PFname(1:4))
    case{'logi'}
        spread = 2.*log((1 - delta - gamma - lambda)./delta)./beta;
    case{'cumu'}
        spread = 2.*sqrt(2).*erfinv(1-(2.*delta./(1 - gamma - lambda)))./beta;
    case{'weib'}
        spread = alpha.*((-log(1 - (1 - lambda - delta - gamma)./(1 - gamma - lambda))).^(1./beta)-(-log(1 - delta./(1 - gamma - lambda))).^(1./beta));
    case{'quic'}
        spread = alpha.*((-log2(1 - (1 - lambda - delta - gamma)./(1 - gamma - lambda))).^(1./beta)-(-log2(1 - delta./(1 - gamma - lambda))).^(1./beta));
    case{'gumb'}
        spread = log10(log(delta./(1 - gamma - lambda))./log(1 - delta./(1 - gamma - lambda)))./beta;
    case{'logq'}
        spread = log10(log(delta./(1 - gamma - lambda))./log(1 - delta./(1 - gamma - lambda)))./beta;
    case{'hype'}
        spread = 2.*log(tan((pi./2).*(1 - lambda - delta - gamma)./(1 - gamma - lambda))./tan((pi./2).*(delta)./(1 - gamma - lambda)))./(pi.*beta);
    otherwise
        warning('PALAMEDES:invalidOption','The function %s is not supported',PFname);
        spread = [];
end