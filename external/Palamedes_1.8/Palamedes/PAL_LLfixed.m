%
% PAL_LLfixed     Log Likelihood associated with 0 df binomial model.
%
% Syntax: LL = PAL_LLfixed(NumPos, OutOfNum, pPos)
%
% Internal Function
%
% Introduced: Palamedes version 1.7.0 (NP)

function [LL] = PAL_LLfixed(NumPos, OutOfNum, pPos)

LL = sum(PAL_nansum(NumPos.*log(pPos)))+sum(PAL_nansum((OutOfNum-NumPos).*log(1 - pPos)));