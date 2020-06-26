function rocArea = computeROCArea(signalMean,noiseMean,commonSd,criteria)
% rocArea = computeROCArea(signalMean,noiseMean,commonSd,criteria)
%
% Compute area under ROC curve given signal mean, noise mean, the common
% standard deviation, and a set of criteria.  Uses the equal-variance
% normal assumption.

% Compute ROC curve
for i = 1:length(criteria)
    [pHitAnalytic(i),pFaAnalytic(i)] = analyticpHitpFa(signalMean,noiseMean,commonSd,criteria(i));
end

% Integrate numerically to get area.  The negative
% sign is because the way the computation goes, the hit rates
% decrease with increasing criteria.
rocArea = -trapz([1 pFaAnalytic 0],[1 pHitAnalytic 0]);

end
