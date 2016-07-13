%
%PAL_SDT_Summ_MultiplePFML_Fit   Fits the parameters of a set of summation 
% psychometric curves using a Maximum Likelihood criterion, based on
% either a probability or additive summation model
%
% Syntax: 
% [gParamsOut, pParamsOut, negLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,gParams,pParams,NumPos,OutOfNum,SummFunc,M,Q)
%
%
%Input:
%
%   'StimLevels': matrix of stimulus levels
%       
%   'gParams' vector of guesses for stimulus level scaling
%       factors gA, gB gC etc. 
%
%   'pParams' vector of gueses for corresponding transducer exponents 
%       pA, pB, pC etc.
%
%   'NumPos' matrix of number of positive responses corresponding to
%    each stimulus level
%   
%   'OutOfNum': number of trials per stimulus level
%
%   'SummFunc': either @PAL_SDT_PS_uneqSLtoPC or @PAL_SDT_AS_uneqSLtoPC;
%
%   'M': number of alternatives/intervals in forced-choice task
%
%   'Q': number of monitored channels
%
%Output:
%
%   'gParamsOut'  Parameter estimates of g, the gain on the stimulus 
%
%   'pParamsOut'  Parameter estimates of p, the transducer exponent 
%
%   'negLL': Negative Log likelihood associated with the fit.
%
%   'exitflag': 1 indicates a succesful fit, 0 indicates fit did not
%       converge (trying again using new initial guesses might help).
%
%   'output': message containing some information concerning fitting
%       process.
%
%Example:
%
% This example uses the routine to simultaneously fit three psychometric 
% functions, one for stimulus A, one for stimulus B and one for stimulus 
% A+B, under the Fixed Attention Window scenario, i.e. when the observer is
% monitoring both A and B channels.  The data is fitted with a probability
% summation model for the situation in which stimulus levels for A and B 
% are not the same
%
% StimLevels(1,1,:)=[1 2 3 4 5 6 7 8]; % for stim A alone
% StimLevels(1,2,:)=[0 0 0 0 0 0 0 0]; % channel B for stim A alone
% StimLevels(2,1,:)=[0 0 0 0 0 0 0 0]; % channel A for stim B alone
% StimLevels(2,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B alone
% StimLevels(3,1,:)=[1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2]; %Stim A in Stim A+B
% StimLevels(3,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B in Stim A+B
% 
% NumPos(1,:)=[44 49 53 60 64 72 73 78];
% NumPos(2,:)=[43 46 56 58 68 70 74 79]; 
% NumPos(3,:)=[44 57 59 68 76 76 79 79];
% 
% OutOfNum(1,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(2,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(3,:) = [80 80 80 80 80 80 80 80];
% 
% SummFunc=@PAL_SDT_PS_uneqSLtoPC;
%
% M=2;  % 2-AFC
% Q=2;  % two monitored channels
%
% gParams=[0.3 0.3]; %Initial guesses for gA and gB
% pParams=[1.5 1.5]; %Initial guesses for pA and pB
%
%Fit data:
%
%   [gParamsOut pParamsOut, negLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,gParams,pParams,NumPos,OutOfNum,SummFunc,M,Q)
%
%
%returns:
% 
% gParamsOut =
% 
%   0.2898    0.2508    
%
% pParamsOut = 
%
%   1.2499    1.4779
% 
% PSnegLL =
% 
%   829.1145
% 
% 
% exitflag =
% 
%      1
% 
% 
% output = 
% 
%        message: 'Search converged successfully. TolX = 1.000000e-06. TolFun = 1.000000e-06.'
%     iterations: 145
%      funcCount: 254
%      algorithm: 'Nelder-Mead simplex direct search'
% 

%Introduced: Palamedes version 1.8.0 (FK & NP)


function [gParamsOut, pParamsOut, negLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,gParams,pParams,NumPos,OutOfNum,SummFunc,M,Q)


options = PAL_minimize('options');   %PAL_minimize search options        

paramsIn=[gParams pParams];


[paramsOut, negLL, exitflag, output] = PAL_minimize(@PAL_SDT_Summ_MultiplePFML_negLL,paramsIn,options,StimLevels,NumPos,OutOfNum,SummFunc,M,Q); 


num=length(paramsOut);

gParamsOut=paramsOut(1:num/2);
pParamsOut=paramsOut(num/2+1:num);




