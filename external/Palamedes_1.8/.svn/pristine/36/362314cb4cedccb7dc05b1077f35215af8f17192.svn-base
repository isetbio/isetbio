%
%PAL_MLDS_negLL (negative) Log Likelihood for MLDS fit.
%
%syntax: [negLL] = PAL_MLDS_negLL(FreeParams, stim, NumGreater, OutOfNum)
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)

function [negLL] = PAL_MLDS_negLL(FreeParams, stim, NumGreater, OutOfNum)

PsiValues = [0 FreeParams(1:length(FreeParams)-1) 1];

if size(stim,2) == 4
    D = (PsiValues(stim(:,2))-PsiValues(stim(:,1)))-(PsiValues(stim(:,4))-PsiValues(stim(:,3)));
end
if size(stim,2) == 3
    D = (PsiValues(stim(:,2))-PsiValues(stim(:,1)))-(PsiValues(stim(:,3))-PsiValues(stim(:,2)));
end
if size(stim,2) == 2
    D = (PsiValues(stim(:,1))-PsiValues(stim(:,2)));
end
    
Z_D = D./FreeParams(length(FreeParams));
pFirst = .5 + .5*(1-erfc(Z_D./sqrt(2)));

negLL = -sum(NumGreater(NumGreater > 0).*log(pFirst(NumGreater > 0)))-sum((OutOfNum(NumGreater < OutOfNum)-NumGreater(NumGreater < OutOfNum)).*log(1-pFirst(NumGreater < OutOfNum)));