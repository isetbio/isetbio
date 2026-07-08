% [dprime,critNorm] = computeDPrimeCritNorm(pHit,pFa)
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