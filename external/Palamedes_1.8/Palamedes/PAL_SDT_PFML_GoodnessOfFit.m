%
%PAL_SDT_PFML_GoodnessOfFit Determines goodness-of-fit of the parameters 
%   of a signal detection theory (SDT) model fitted to psychometric 
%   function (PF) data
%
%   Note that depending on the number of simulations and the particular SDT 
%   model the function could take some time to execute
%
% Syntax: 
%   [Dev,pDev,DevSim,converged] = PAL_SDT_PFML_GoodnessOfFit(StimLevels,...
%       params,NumPosData,OutOfNum,SDTfunc,M,B)
%
%Input:
%
%   'StimLevels': vector of stimulus levels
%   
%   'params': 2-valued vector of parameters. These are the stimulus level 
%       scaling factor g and the transducer exponent p.
%
%   'OutOfNum': number of trials per stimulus level
%
%   'SDTfunc': examples @PAL_SDT_2AFC_DPtoPC, 
%       @PAL_SDT_3AFCoddity_IndMod_DPtoPC
%
%   'M': number of alternatives in forcd-choice task.  If SDTfunc does not
%       take M as an argument M should be empty, i.e. 'M=[]'
%
%   'B': number of simulations
%
%Output:
%
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric functionx to fit of saturated model)
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater than Deviance value of data. The greater the value of pDev,
%       the better the fit. Typically fit can be rejected if pDev<0.05
%
%   'DevSim': vector containing all B simulated Deviance values.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%
%Example:
%
%   StimLevels=[1 2 3 4 5 6];  
%   params=[0.3222 1.8320];
%   NumPosData=[51 60 79 85 96 99];
%   OutOfNum=[100 100 100 100 100 100];
%   SDTfunc=@PAL_SDT_2AFC_DPtoPC;
%   M=[];
%   B=400;
%
%Fit data:
%
%   [Dev, pDev] = PAL_SDT_PFML_GoodnessOfFit(StimLevels,params,NumPosData,OutOfNum,SDTfunc,M,B)
%
%returns something like:
%
% Dev =
% 
%     1.9372
% 
% pDev =
% 
%     0.8000
%
%Introduced: Palamedes version 1.8.0 (FK & NP)


function [Dev, pDev, DevSim, converged] = PAL_SDT_PFML_GoodnessOfFit(StimLevels,params,NumPosData,OutOfNum,SDTfunc,M,B)

converged = false(B,1);
DevSim = zeros(B,1);
pDev = [];

negLLCon = PAL_SDT_PFML_negLL(params,StimLevels,NumPosData,OutOfNum,SDTfunc,M);
negLLAug = PAL_PFML_negLLNonParametric(NumPosData,OutOfNum);

Dev = 2*(negLLCon-negLLAug);

for b = 1:B
    
    NumPosSim = PAL_SDT_PF_SimulateObserverParametric(StimLevels,params,OutOfNum,SDTfunc,M);
    
    [trash, negLLConSim, converged(b)] = PAL_SDT_PFML_Fit(StimLevels,NumPosSim,OutOfNum,params,SDTfunc,M);
    
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