%
%PAL_SDT_1AFC_ROCML_negLLNonParametric     (negative) Log Likelihood
% associated with saturated model for a 1AFC ROC curve
%
%
%Internal Function
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [negLL, numParams] = PAL_SDT_ROCML_negLLNonParametric(cumNumHF, OutOfNum) 

pHF=cumNumHF./OutOfNum;

negLL = -sum(PAL_nansum(cumNumHF.*log(pHF))+PAL_nansum((OutOfNum-cumNumHF).*log(1 - pHF)));
numParams = sum(sum(OutOfNum~=0));