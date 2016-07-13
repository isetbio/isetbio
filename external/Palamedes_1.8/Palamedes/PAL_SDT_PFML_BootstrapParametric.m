%
%PAL_SDT_PFML_BootstrapParametric    Perform parametric bootstrap to
%   determine standard errors on parameters of an SDT (signal detection
%   theory) model fitted to PF (psychometric function) data
%
%   Note that depending on the number of simulations and the particular SDT 
%   model the function could take some time to execute
%
% Syntax: 
% [gSE pSE] = PAL_SDT_PFML_BootstrapParametric(StimLevels,params,...
%   OutOfNum,SDTfunc,M,B)
%
%
%Input:
%
%   'StimLevels': vector of stimulus levels
%   
%   'params': 2-valued vector of parameters. These are the stimulus level 
%       scaling factor g and the transducer exponent p.
%
%   'paramsFree': 3-valued vector with 1s for free and 0s for fixed 
%       parameters;
%
%   'OutOfNum': number of trials per stimulus level
%
%   'SDTfunc': examples @PAL_SDT_2AFC_DPtoPC, 
%       @PAL_SDT_3AFCoddity_IndMod_DPtoPC;
%
%   'M': number of alternatives in forcd-choice task.  If SDTfunc does not
%       take M as an argument M should be empty, i.e. 'M=[]'
%
%   'B': number of bootstrap simulations
%
%Output:
%
%   'gSE': scalar standard error of estimates of stimulus scaling factor g
%
%   'pSE': scalar standard error of estimates of transducer exponent p
%
%
%Example:
%
%   StimLevels=[1 2 3 4 5 6];  
%   params=[0.3222 1.8320];
%   OutOfNum=[100 100 100 100 100 100];
%   SDTfunc=@PAL_SDT_2AFC_DPtoPC;
%   M=[];
%   B=400;
%
%Fit data:
%
%   [gSE pSE] = PAL_SDT_PFML_BootstrapParametric(StimLevels,params,...
%       OutOfNum,SDTfunc,M,B)
%
%
%returns something like:
% 
% gSE =
% 
%     0.0216
% 
% 
% pSE =
% 
%     0.3125
% 
%Introduced: Palamedes version 1.8.0 (FK & NP)


function [gSE, pSE, alphaSE] = PAL_SDT_PFML_BootstrapParametric(StimLevels,params,OutOfNum,SDTfunc,M,B)

converged = false(B,1);
g=zeros(B,1);
p=zeros(B,1);

for b = 1:B

    %Simulate experiment
    NumPos = PAL_SDT_PF_SimulateObserverParametric(StimLevels,params,OutOfNum,SDTfunc,M);
    
    %fit data
    [paramsOut, trash, converged(b)] = PAL_SDT_PFML_Fit(StimLevels,NumPos,OutOfNum,params,SDTfunc,M);
    
    g(b)=paramsOut(1);
    p(b)=paramsOut(2);
  
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