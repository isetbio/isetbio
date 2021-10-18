function [fractionCorrect,dPrime,decisionMean,decisionVar] = analyticPoissonIdealObserver(alphaMeanResponses,betaMeanResponses)
%analyticPoissonIdealObserver]  Get the ideal observer fraction correct in a TAFC task
%    [fractionCorrect,dPrime,decisionMean,decisionVar] = analyticPoissonIdealObserver(alphaMeanResponses,betaMeanResponses)
%
%    Get the fraction correct in a TAFC task using the formula for a
%    Poisson ideal observer developed by Geisler, 1984, JOSA A, 1, pp. 775
%    ff.
%
%    Also returns dPrime.  Indeed, the Geisler formula returns dPrime and
%    this routine then converts to TAFC percent correct by numerically
%    integrating the area under the ROC curve for that dPrime, assuming
%    normal distributions with equal variance in the conversion.
%
%    The return variables decisionMean and decisionVar are the mean and
%    variance of the normal approximation to the decision
%    variable.
%
%    The two input arguments should all be vectors of the same length, and
%    give the mean of the Poisson distributed responses for the two
%    stimulus types.
%
% See also dPrimeToTAFCFractionCorrect

% Handle special case when the mean responses for the two classes are
% identical, which runs into numerical problems if we compute it using the
% code below
if (all(alphaMeanResponses == betaMeanResponses))
    fractionCorrect = 0.5;
    dPrime = 0;
    decisionMean = 0;
    decisionVar = 1;
    return;
end

% Sometimes the mosaic extends beyond the stimulus, giving 0 mean responses
% and causing trouble in the formula.  In general, we need only consider
% locations where the alpha and beta means differ, so let's handle that.
index = find(alphaMeanResponses ~= betaMeanResponses);

% This comes from the appendix of the paper
decisionMean = sum( (betaMeanResponses(index)-alphaMeanResponses(index)).*log(betaMeanResponses(index)./alphaMeanResponses(index)) );
decisionVar = 0.5*sum( (betaMeanResponses(index)+alphaMeanResponses(index)).*(log(betaMeanResponses(index)./alphaMeanResponses(index)).^2) );
dPrime = decisionMean / sqrt(decisionVar);

% Call into our dPrime to TAFC conversion function
fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime);

if (isnan(dPrime) || isnan(fractionCorrect))
    error('Should not get NaN here.  Probably there are some zero mean response cones');
end



