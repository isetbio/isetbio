%
%PAL_SDT_PFML_Fit   Fits the parameters of an SDT (signal dtection theory)
%   psychometric function (PF) using a Maximum Likelihood criterion.
%
% Syntax: 
% [paramsOut,negLL,exitflag,output] = PAL_SDT_Summ_PFML_Fit(StimLevels,...
%   NumPos,OutOfNum,params,SDTfunc,M)
%
%
%Input:
%
%   'StimLevels': vector of stimulus levels
%       
%   'NumPos': vector of number of positive responses corresponding to
%    each stimulus level
%   
%   'OutOfNum': number of trials per stimulus level
%
%   'params': 2-valued vector of parameters.  These are the guesses 
%       for the free parameters stimulus level scaling factor g and 
%       transducer exponent p.  
%
%   'SDTfunc': examples @PAL_SDT_2AFC_DPtoPC, 
%       @PAL_SDT_3AFCoddity_IndMod_DPtoPC;
%
%   'M': number of alternatives in forcd-choice task.  If SDTfunc does not
%       take M as an argument, M should be empty, i.e. 'M=[]'
%
%Output:
%
%   'paramsOut'  Parameter estimates of g and p 
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
%   StimLevels=[1 2 3 4 5 6];  
%   params=[0.25 1.25];
%   NumPos=[51 60 79 85 96 99];
%   OutOfNum=[100 100 100 100 100 100];
%   SDTfunc=@PAL_SDT_2AFC_DPtoPC;
%   M=[];
%
%Fit data:
%
%   [paramsOut, negLL, exitflag, output] = ...
%       PAL_SDT_PFML_Fit(StimLevels,NumPos,OutOfNum,params,SDTfunc,M)
%
%returns:
% 
% paramsOut =
% 
%     0.3222    1.8380
% 
% negLL =
% 
%   253.6256
% 
% exitflag =
% 
%      1
%
% output = 
% 
%      message: 'Search converged successfully. TolX = 1.000000e-06. 
%           TolFun = 1.000000e-06.'
%      iterations: 55
%      funcCount: 106
%      algorithm: 'Nelder-Mead simplex direct search'
%
%Introduced: Palamedes version 1.8.0 (FK & NP)


function [paramsOut, negLL, exitflag, output] = PAL_SDT_PFML_Fit(StimLevels,NumPos,OutOfNum,params,SDTfunc,M)

options = []; %default

[paramsOut, negLL, exitflag, output] = PAL_minimize(@PAL_SDT_PFML_negLL,params,options,StimLevels,NumPos,OutOfNum,SDTfunc,M);