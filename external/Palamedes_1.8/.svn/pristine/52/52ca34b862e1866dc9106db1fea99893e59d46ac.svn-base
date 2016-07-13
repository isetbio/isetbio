%
% PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a 2AFC (two-alternative-forced-choice) match-to-sample
% task, assuming a Differencing Observer model and an unbiased 
% observer
%
% Syntax: [dP]=PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP(pC);
% 
% returns a matrix of d' ('dP') for an input matrix of proportion correct 
% ('pC'), defined in the range 0<p<1 (p=proportion)
%
% Example:
%
% [dP]=PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP([.5 .6 .7 .8 .9])
%
% returns:
% 
% dP =
% 
%    -0.0000    1.1152    1.7153    2.3549    3.2631
%
% The example input is an N=5 vector of proportion correct and the output 
% the corresponding d's
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)

function [dP]=PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP(pC)

func=@PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC;

dP = zeros(1,length(pC));

for r=1:length(pC)
    dP(r)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pC(r),func,[]);
end