%
% PAL_SDT_2AFCmatchSample_IndMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a 2AFC (two-alternative-forced-choice) match-to-sample
% task, assuming an Independent Observer model and an unbiased 
% observer
%
% Syntax: [dP]=PAL_SDT_2AFCmatchSample_IndMod_PCtoDP(pC);
% 
% returns a scalar, vector or matrix of d' ('dP') for a scalar, vector or
% matrix of proportion correct p ('pC'), defined in the range 0<p<1 
%
% Example:
%
% [dP]=PAL_SDT_2AFCmatchSample_IndMod_PCtoDP([.5 .6 .7 .8 .9])
%
% returns:
% 
% dP =
% 
%    -0.0000    1.0020    1.5287    2.0726    2.8073
%
% The example input is an N=5 vector of proportion correct and the output 
% the corresponding d's
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)

function [dP]=PAL_SDT_2AFCmatchSample_IndMod_PCtoDP(pC)

func=@PAL_SDT_2AFCmatchSample_IndMod_DPtoPC;

dP = zeros(1,length(pC));

for r=1:length(pC)    
    dP(r)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pC(r),func,[]);
end