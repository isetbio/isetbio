function dPrime = TAFCFractionCorrectToDPrime(fractionCorrect)
% dPrime = TAFCFractionCorrectToDPrime(fractionCorrect)
%
% Convert TAFC fraction correct to an equal-variance normal distribution
% d-prime.
%
% See also: dPrimeToTAFCFractionCorrect, computeROCArea, analyticPHitPpFA, computeDPrimCritNorm

% History:
%    05/26/2020  dhb  Wrote it.

%% Examples:
%{
    fractionCorrect = dPrimeToTAFCFractionCorrect(1);
    dPrime = TAFCFractionCorrectToDPrime(fractionCorrect);
    if (max(abs(dPrime-1)) > 1e-4)
        error('d-prime functions do not self invert');
    end
%}

dPrime = sqrt(2)*norminv(fractionCorrect);
    
end





