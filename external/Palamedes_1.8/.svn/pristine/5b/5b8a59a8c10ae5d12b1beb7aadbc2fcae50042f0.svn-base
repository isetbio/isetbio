%
% PAL_SDT_3AFCoddity_IndMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a 3-AFC (3-alternative-forced-choice) oddity
% task, assuming an Independent Observation model and an unbiased observer
%
% Syntax: [dP] = PAL_SDT_MAFCoddity_IndMod_PCtoDP(pC,{optional arguments});
% 
% returns a scalar, vector or matrix of d' ('dP') for an input scalar, 
% vector or matrix of proportion correct p ('pC'), defined in the ranges 
% 0<p<1. 
%
% Note that unlike PAL_SDT_MFCoddity_IndMod_PCtoDP, which is the  
% M-AFC version of the same task but uses an inverse routine that uses
% Monte Carlo simulation, the inverse of this routine is deterministic and 
% hence this routine runs much faster
%
% PAL_SDT_3AFCoddity_IndMod_PCtoDP uses Nelder-Mead Simplex method as 
%   implemented in PAL_minimize. The default search options may be changed 
%   by using the following syntax:
%
%   [dP] = PAL_SDT_MAFCoddity_IndMod_PCtoDP(pC,'searchOptions', options) 
%
%   where 'options' is a structure created using: options = 
%   PAL_minimize('options'). For more information type 
%   PAL_minimize('options','help').
%
%
% Example:
%
% [dP] = PAL_SDT_3AFCoddity_IndMod_PCtoDP([.5 .6 .7 .8 .9])
%
% returns:
% 
% dP =
% 
%     1.2451    1.6686    2.1012    2.6060    3.3245
%
% The example input arguments are an N=5 vector of proportion correct and 
% the output is an N=5 vector of d's
%
% Introduced: Palamedes version 1.6.0 (FK & NP)
% Modified: Palamedes version 1.6.3 (see History.m)

function [dP] = PAL_SDT_3AFCoddity_IndMod_PCtoDP(pC, varargin)

options = []; %default

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end

[rows, cols] = size(pC);

func = @PAL_SDT_3AFCoddity_IndMod_DPtoPC;

dP = zeros(rows, cols);

for r = 1:rows
    for c = 1:cols
        dP(r,c) = PAL_minimize(@PAL_sqDistanceYfuncX,1,options,pC(r,c),func,[]);
        if pC(r,c) == 1.0
            dP(r,c) = 1/0;
        end
    end
end