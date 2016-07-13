%
%PAL_SDT_SLtoPC   Converts stimulus intensities to proportion correct 
%   assuming parameters of an SDT (signal detection theory) model
%
% Syntax: 
%
% [PC] = PAL_SDT_SLtoPC(x,g,p,SDTfunc,M)
%
%Input:
%
%   'x': scalar, vector or matrix of stimulus levels
%
%   'g': stimulus level scaling factor 
%
%   'p': transducer exponent p  
%
%   'SDTfunc': SDT routine that converts d-prime to prop. correct,
%       examples @PAL_SDT_2AFC_DPtoPC, @PAL_SDT_3AFCoddity_IndMod_DPtoPC;
%
%   'M': number of alternatives in forced-choice task. If SDTfunc does
%       not take M as an argument, M should be set to empty, i.e. 'M=[]'
%
%Output:
%
%   'PC':  scalar, vector or matrix of prop. correct
%
%Example:
%
%   x=[1 2 3 4 5];  
%   g=0.5;
%   p=1.5;
%   SDTfunc=@PAL_SDT_2AFC_DPtoPC;
%   M=[];
%
%Calculate proportion correct:
%
%   [PC] = PAL_SDT_SLtoPC(x,g,p,SDTfunc,M)
%
%returns:
% 
% PC =
% 
%   0.5987    0.7602    0.9030    0.9772    0.9974  
% 
%Introduced: Palamedes version 1.8.0 (FK)


function [PC] = PAL_SDT_SLtoPC(x,g,p,SDTfunc,M)

DP=(g.*x).^p;

%obtain proportion correct from d'prime
if isempty(M)
    PC = SDTfunc(DP); 
else
    PC = SDTfunc(DP,M);
end





