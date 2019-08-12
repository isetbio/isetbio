function fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime)
% fractionCorrect = dPrimeToTAFCPercentFraction(dPrime)
%
% Get area under ROC curve for a normal distribution d-prime and from there
% compute fraction correct.

nCriteria = 1000;
lowCriterion = -8;
highCriterion = 8;

fractionCorrect = ComputeROCArea(dPrime,0,1,linspace(lowCriterion,highCriterion,nCriteria));
    
end

%% Function to get area under ROC curve
function rocArea = ComputeROCArea(signalMean,noiseMean,commonSd,criteria)

% Compute ROC curve
for i = 1:length(criteria)
    [pHitAnalytic(i),pFaAnalytic(i)] = AnalyticpHitpFa(signalMean,noiseMean,commonSd,criteria(i));
end

% Integrate numerically to get area.  The negative
% sign is because the way the computation goes, the hit rates
% decrease with increasing criteria.
rocArea = -trapz([1 pFaAnalytic 0],[1 pHitAnalytic 0]);

end


%% [pHit,pFa] = AnalyticpHitpFa(signalMean,noiseMean,commonSd,rightCrit)
% 
% This just finds the area under the signal and noise
% distributions to the right of the criterion to obtain
% hit and false alarm rates
function [pHit,pFa] = AnalyticpHitpFa(signalMean,noiseMean,commonSd,rightCrit)

pHit = 1-normcdf(rightCrit,signalMean,commonSd);
pFa = 1-normcdf(rightCrit,noiseMean,commonSd);

end

%% [dprime,critNorm] = ComputeDPrimeCritNorm(pHit,pFa)
%
% FindFdPrime and criterion from
% hit and fa rates.
%
% This assumes equal variance normal for the noise and
% signal response distributions.
%
% The criterion is returned in normalized units where the
% noise distribution is taken to have mean 0 and the common SD is 1.
%
% Formula from lecture slides.
function [dprime,critNorm] = ComputeDPrimeCritNorm(pHit,pFa)

dprime = norminv(pHit,0,1)-norminv(pFa,0,1);
critNorm = norminv(1-pFa,0,1);

end

