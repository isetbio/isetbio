%
%PAL_SDT_PCtoSL   Converts proportion correct to stimulus intensity 
%   assuming parameters of an SDT (signal detection theory) model
%
% Syntax: 
%
% [x] = PAL_SDT_PCtoSL(PC,g,p,invSDTfunc,M)
%
%Input:
%
%   'PC': scalar, vector or matrix of proportion correct detections 
%
%   'g': stimulus level scaling factor 
%
%   'p': transducer exponent p  
%
%   'invSDTfunc': SDT routine that converts prop. correct to d-prime,
%       examples @PAL_SDT_2AFC_PCtoDP, @PAL_SDT_3AFCoddity_IndMod_PCtoDP;
%
%   'M': number of alternativs in forced-choice task. If invSDTfunc does
%       not take M as an argument, M must be set to empty, i.e. 'M=[]'
%
%Output:
%
%   'x'  scalar, vector or matrix of stimulus intensity
%
%
%Example:
%
%   PC=[0.55 0.6 0.65 0.7 0.9];  
%   g=3;
%   p=1.5;
%   invSDTfunc=@PAL_SDT_2AFC_PCtoDP;
%   M=[];
%
%Calculate stimulus levels:
%
%   [x] = PAL_SDT_PCtoSL(PC,g,p,invSDTfunc,M)
%
%returns:
% 
% x =
% 
%     0.1054    0.1682    0.2224    0.2731    0.4955
%    
%Introduced: Palamedes version 1.8.0 (FK)

function [x] = PAL_SDT_PCtoSL(PC,g,p,invSDTfunc,M)


%obtain d'prime at prop. correct level 
if isempty(M)
    DP = invSDTfunc(PC); 
else
    DP = invSDTfunc(PC,M);
end

x=(1/g).*DP.^(1/p);