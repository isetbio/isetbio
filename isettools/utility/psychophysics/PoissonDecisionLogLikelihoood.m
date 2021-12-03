% Log-likelihood for Poisson response vectors 
%
% Syntax:
%    llDecision = PoissonDecisionLogLikelihoood(response, template)
%
% Description:
%    Compute the part of the Poisson log likelihood that depends
%    on the template, for use in Poisson ideal observer decisions.
%
%    If we have Poisson responses (response) and the Poisson mean
%    template, then the PMF is
%        p = template.^response*exp(-template)/(response!)
%    Take the log and we have
%        log(p) = response*log(template) - template - log(response!)
%
%    Computing log(response!) is slow for large response, but we don't need
%    that if all we're going to do is compare log likelihoods across
%    different templates, for fixed response.  So here we just compute and
%    return the first two terms.
%
%    
%    Poisson observations, given Poisson mean for each.
%    See for example:
%      https://online.stat.psu.edu/stat504/node/28/
%
% Inputs:
%    response  -   Vector of observed responses.
%    template  -   Poisson mean for each entry of response.
%
% Outputs:
%    llDecision -   Log likelihood, without -log(response!) term.
%
% See also: PoissonIdealObserverNAlternativeFC
%

% History:
%    10/17/20  dhb  Added comments.
%    12/03/21  dhb  Rename, pull out, better comments.

function llDecision = PoissonDecisionLogLikelihoood(response, template)
    response = double(response(:));
    template = double(template(:));
    if (length(response) ~= length(template))
        error('Passed response and template must have same number of entries');
    end
    llDecision = sum(response .* log(template) - template);
end
