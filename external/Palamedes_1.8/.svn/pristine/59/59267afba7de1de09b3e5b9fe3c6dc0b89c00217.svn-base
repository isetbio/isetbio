%
%PAL_SDT_Summ_MultiplePFML_BootstrapParametric    Perform parametric 
%   bootstrap to determine standard errors on parameters of SDT 
%   (signal detection theory) probability and additive summation models 
%   simultaneously fitted to multiple PF (psychometric function) 
%   data.  
%
%   Note that depending on the number of simulations and the number of PFs 
%   the routine may take some time to execute.
%   
%   
% Syntax: 
% [gSE, pSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q,B)
%
%
%Input:
%
%   'StimLevels': matrix of stimulus levels
%       
%   'gParams' vector of fitted g (stimulus gain) parameters
%
%   'pParams' vector of fitted p (transducer exponent) parameters
%   
%   'OutOfNum': number of trials per stimulus level
%
%   'SummFunc': @PAL_SDT_PS_uneqSLtoPC, @PAL_SDT_AS_uneqSLtoPC,
%               @PAL_SDT_PS_SLtoPC or @PAL_SDT_AS_SLtoPC;
%
%   'M': number of alternatives/intervals in forced-choice task
%
%   'Q': number of monitored channels
%
%   'B': number of simulations
%
%Output:
%
%   'gSE'  standard error of g parameter estimate
%
%   'pSE'  standard error of p parameter estimates 
%
%
%Example:
%
% This example uses the routine to estimate the errors on p and g 
% when these have been estimated from the simultaneous fit of three 
% psychometric functions, one for stimulus A, one for stimulus B and one 
% for stimulus A+B, under the Fixed Attention Window scenario, i.e. when 
% the observer is monitoring both A and B channels.  The data is fitted 
% with a probability summation model for the situation in which stimulus 
% levels for A and B are not the same
%
% StimLevels(1,1,:)=[1 2 3 4 5 6 7 8]; % for stim A alone
% StimLevels(1,2,:)=[0 0 0 0 0 0 0 0]; % channel B for stim A alone
% StimLevels(2,1,:)=[0 0 0 0 0 0 0 0]; % channel A for stim B alone
% StimLevels(2,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B alone
% StimLevels(3,1,:)=[1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2]; %Stim A in Stim A+B
% StimLevels(3,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B in Stim A+B
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
% B=200; % 200 simulations to obtain standard error
%
% gParams=[0.2998 0.2508]  %fitted g parameters for the two channels
% pParams=[1.2499 1.4779]; %fitted p parameters for the two channels
%
%Obtain bootstrap errors:
%
%  [gSE, pSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q,B)
% 
%
%returns something like:
% 
% gSE =
% 
%     0.0249    0.0194
% 
% 
% pSE =
% 
%     0.1760    0.2141
%  
%Introduced: Palamedes version 1.8.0 (FK & NP)


function [gSE, pSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q,B)


converged = false(B,1);


num=length(gParams);
g=zeros(B,num);
p=zeros(B,num);


for b = 1:B
    
    %Simulate experiment
    NumPos = PAL_SDT_Summ_MultiplePF_SimulateObserverParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q);
    
    [gParamsOut, pParamsOut, trash, converged(b)] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,gParams,pParams,NumPos,OutOfNum,SummFunc,M,Q);
    
    g(b,:)=gParamsOut(:);
    p(b,:)=pParamsOut(:);
    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end
    
end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end


[trash, gSE] = PAL_MeanSDSSandSE(g);
[trash, pSE] = PAL_MeanSDSSandSE(p);
