%
% PAL_PFML_GoodnessOfFitZeroDF     Determine Goodness-of-Fit of model of
%   binomial data where the model consists only of fixed parameter values
%   (i.e., the model has zero degrees of freedom).
%
% Syntax:   [Dev pDev DevSim converged] = ...
%           PAL_GoodnessOfFitZeroDF(NumPos, OutOfNum, pPos, B);
%
% Input:
%   'NumPos': array containing frequency of positive response (e.g., 
%       'correct') for each trial type used.
%
%   'OutOfNum': array containing for each of the entries of 'NumPos' the 
%       number of trials used.
%
%   'pPos': array containing for each of the entries of 'NumPos' the 
%       probability of a positive response according to the model tested.
%       May also be scalar if all equal.
%
%   'B': number of Monte Carlo simulations to perform. May be set to 0 in
%       case one prefers to use asymptotic chi-square distribution to
%       determine p-value using e.g., Matlab function chi2cdf in the
%       (optional) Matlab 'statistics' toolbox.
%
%Output: 
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       model to fit of saturated model). Asymptotically distributed as
%       chi-square with degrees of freedom equal to number of conditions.
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater or equal to Deviance value of data. The greater the value 
%       of pDev, the better the fit.
%
%   'DevSim': vector containing all B simulated Deviance values.
%
%   'converged': Will always be logical(true). Included only for analogy to
%       other Goodness-Of-Fit routines.
%   
% Example 1: 100 flips of coin result in 61 heads, 36 tails. Fair?
%
%   [Dev pDev] = PAL_GoodnessOfFitZeroDF(61,100,.5,1000)
%
% might return:
%   Dev = 
%       4.8798
%   pDev =
%       0.0349
%
% Example 2: Coin 1 flipped 40 times, 25 heads, coin 2 flipped 60 times, 34
%   heads. Both coins fair?
%
%   [Dev pDev] = PAL_GoodnessOfFitZeroDF([25 34],[40 60],.5,1000)
%
% might return:
%   Dev = 
%       3.5966
%   pDev =
%       0.1620
%
% Introduced: Palamedes version 1.7.0 (NP)

function [Dev, pDev, DevSim, converged] = PAL_GoodnessOfFitZeroDF(NumPos, OutOfNum, pPos, B)

converged = logical(true);
DevSim = zeros(B,1);

LLsat = PAL_PFML_LLsaturated(NumPos, OutOfNum);
LLfixed = PAL_LLfixed(NumPos,OutOfNum,pPos);
Dev = 2*(LLsat-LLfixed);

DevSim = zeros(1,B);
trash = ones(size(NumPos));
for rep = 1:B
    NumPos = PAL_PF_SimulateObserverNonParametric(trash, pPos.*OutOfNum, OutOfNum);
    LLsat = PAL_PFML_LLsaturated(NumPos, OutOfNum);
    LLfixed = PAL_LLfixed(NumPos,OutOfNum,pPos);
    DevSim(rep) = 2*(LLsat-LLfixed);
end
pDev = length(DevSim(DevSim>=Dev))/B;