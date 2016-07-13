%
%PAL_PFML_negLLNonParametric     (negative) Log Likelihood associated with 
%   saturated model.
%
%Syntax: negLL = PAL_negLLNonParametric(NumPos, OutOfNum)
%
%Requires trials to have been grouped (e.g., using PAL_PFML_GroupTrialsByX)
%
%Internal Function
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.1.0, 1.6.0, 1.6.3 (see History.m)

function [negLL, numParams] = PAL_PFML_negLLNonParametric(NumPos, OutOfNum)

pcorrect = NumPos./OutOfNum;

negLL = -sum(PAL_nansum(NumPos.*log(pcorrect)))-sum(PAL_nansum((OutOfNum-NumPos).*log(1 - pcorrect)));
numParams = sum(sum(OutOfNum~=0));