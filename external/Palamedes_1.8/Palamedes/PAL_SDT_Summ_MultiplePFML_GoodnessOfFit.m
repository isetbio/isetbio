%
%PAL_SDT_Summ_MultiplePFML_GoodnessOfFit  Obtains Goodness-of-fit of 
%   parameters of SDT (signal detection theory) probability and 
%   additive summation models obtained from simultaneously fitted multiple 
%   PF (psychometric function) data.  
%
%   Note that depending on the number of simulations and the number of PFs 
%   the routine may take some time to execute.
%
% Syntax: 
% [Dev, pDev, DevSim, converged] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,gParams,pParams,NumPosData,OutOfNum,SummFunc,M,Q,B)
%
%
%Input:
%
%   'StimLevels': matrix of stimulus levels
%       
%   'gParams': vector of fitted g (stimulus gain) parameters
%
%   'pParams': vector of fitted p (transducer exponent) parameters
%
%   'NumPosData': matrix of number of correct responses
%   
%   'OutOfNum': matrix of number of trials per stimulus level
%
%   'SummFunc': @PAL_SDT_PS_uneqSLtoPC, @PAL_SDT_AS_uneqSLtoPC,
%               @PAL_SDT_PS_SLtoPC or @PAL_SDT_AS_SLtoPC;
%
%   'M': scalar number of alternatives/intervals in forced-choice task
%
%   'Q': scalar number of monitored channels
%
%   'B': scalar number of simulations
%
%Output:
%
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric functionx to fit of saturated model)
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater than Deviance value of data. The greater the value of pDev,
%       the better the fit.
%
%   'DevSim': vector containing all B simulated Deviance values.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%
%
%Example:
%
% This example uses the routine to estimate the goodness-of-fit
% of a probability summation model applied to three simultaneously fitted 
% psychometric functions, one for stimulus A, one for stimulus B and one 
% for stimulus A+B, under the Fixed Attention Window scenario, i.e. when 
% the observer is monitoring both A and B channels.  The summation model 
% assumes that the stimulus levels for A and B are not the same
%
% StimLevels(1,1,:)=[1 2 3 4 5 6 7 8]; % for stim A alone
% StimLevels(1,2,:)=[0 0 0 0 0 0 0 0]; % channel B for stim A alone
% StimLevels(2,1,:)=[0 0 0 0 0 0 0 0]; % channel A for stim B alone
% StimLevels(2,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B alone
% StimLevels(3,1,:)=[1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2]; %Stim A in Stim A+B
% StimLevels(3,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B in Stim A+B
%
% NumPosData(1,:)=[44 49 53 60 64 72 73 78];
% NumPosData(2,:)=[43 46 56 58 68 70 74 79]; 
% NumPosData(3,:)=[44 57 59 68 76 76 79 79];
% 
% OutOfNum(1,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(2,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(3,:) = [80 80 80 80 80 80 80 80];
% 
% SummFunc=@PAL_SDT_PS_uneqSLtoPC; % probability summation model for 
%       unequal stimulus levels
%
% M=2;  % 2-AFC
% Q=2;  % two monitored channels
% B=200; % 200 simulations to obtain goodness-of-fit
%
% gParams=[0.2998 0.2508]; fitted g parameters for two channels from PS model
% pParams=[1.2499 1.4779]; fitted p parameters for two channels from PS model
%
%Obtain goodness-of-fit:
%
%  [Dev, pDev] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,gParams,pParams,NumPosData,OutOfNum,SummFunc,M,Q,B)
%
%returns something like:
% 
% Dev =
% 
%     7.2422
% 
% 
% pDev =
% 
%      0.9933
% 
%Introduced: Palamedes version 1.8.0 (FK & NP)


function [Dev, pDev, DevSim, converged] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,gParams,pParams,NumPosData,OutOfNum,SummFunc,M,Q,B)


converged = false(B,1);
DevSim = zeros(B,1);
pDev = [];


params=[gParams pParams];


negLLCon = PAL_SDT_Summ_MultiplePFML_negLL(params,StimLevels,NumPosData,OutOfNum,SummFunc,M,Q);
negLLAug = PAL_PFML_negLLNonParametric(NumPosData,OutOfNum);

Dev = 2*(negLLCon-negLLAug);          

for b = 1:B
    
    NumPosSim = PAL_SDT_Summ_MultiplePF_SimulateObserverParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q);
    
    [trash, trash, negLLConSim, converged(b)] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,gParams,pParams,NumPosSim,OutOfNum,SummFunc,M,Q);
    negLLAugSim = PAL_PFML_negLLNonParametric(NumPosSim,OutOfNum);
    
    DevSim(b) = 2*(negLLConSim-negLLAugSim);
    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end
    
    
end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

if B > 0
    pDev = length(DevSim(DevSim>Dev))/B;
end