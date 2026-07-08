function [pHit,pFa] = analyticpHitpFa(signalMean,noiseMean,commonSd,rightCrit)
% [pHit,pFa] = analyticpHitpFa(signalMean,noiseMean,commonSd,rightCrit)
% 
% This just finds the area under the signal and noise
% distributions to the right of the criterion to obtain
% hit and false alarm rates.  Uses the equal-variance normal assumption.

pHit = 1-normcdf(rightCrit,signalMean,commonSd);
pFa = 1-normcdf(rightCrit,noiseMean,commonSd);

end

