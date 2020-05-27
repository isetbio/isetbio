function fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime)
% fractionCorrect = dPrimeToTAFCPercentFraction(dPrime)
%
% Get area under ROC curve for an equal-variance normal distribution
% d-prime and from there compute fraction correct.
%
% This was originally written as numerical integration of ROC curve.  Now
% has analytic calculation inserted, and a check that the two are close.
% At some point, might want to change over to the analytic version.
%
% See also: computeROCArea, analyticPHitPpFA, computeDPrimCritNorm

% History:
%    05/26/2020  dhb  Added analytic calculation

%% Examples:
%{
    dPrimeToTAFCFractionCorrect(0.25)
    dPrimeToTAFCFractionCorrect(0.5)
    dPrimeToTAFCFractionCorrect(0.1)
    dPrimeToTAFCFractionCorrect(2)
    dPrimeToTAFCFractionCorrect(3)
%}

nCriteria = 1000;
lowCriterion = -8;
highCriterion = 8;

fractionCorrect = computeROCArea(dPrime,0,1,linspace(lowCriterion,highCriterion,nCriteria));
fractionCorrect1 = normcdf(dPrime/sqrt(2));
if (max(abs(fractionCorrect-fractionCorrect1)) > 1e-3)
    error('Numerical and analytic calculations do not match');
end
    
end





