%
% PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP converts proportion hits and 
% proportion false alarms into d'(d-prime) and criterion C for a 2AFC 
% (two-alternative-forced-choice) match-to-sample task under the
% Independent Observation model
%
% Syntax: [dP C]=PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP(pHF);
% 
% returns a scalar or N-length vector of d' ('dP') and criterion C ('C') 
% for an Nx2 input matrix of N proportion hits and proportion false alarms 
% ('pHF') defined in the range 0<p<1 (p=proportion)
%
% Example:
%
% [dP C]=PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP([.6 .1; .8 .1; .9 .6])
%
% returns:
% 
% dP =
% 
%     1.9481
%     2.4390
%     1.5103
% 
% 
% C =
% 
%     0.5141
%     0.2200
%    -0.7674
%
% The example input argument is a 3 x 2 matrix with each row (demarcated
% by a semi-colon) consisting of a proportion of hits and proportion of 
% false alarms.  The columns in the output are the resulting N=3 
% vectors of dP and C
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)

function [dP, C]=PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP(pHF)

[rows, cols]=size(pHF);

zH=PAL_PtoZ(pHF(:,1));
zF=PAL_PtoZ(pHF(:,2));

C=-0.5.*(zH+zF);

zDiff=(zH-zF)./2;
PCmax=PAL_ZtoP(zDiff);

func=@PAL_SDT_2AFCmatchSample_IndMod_DPtoPC;

dP = zeros(rows,1);

for r=1:rows
    dP(r)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],PCmax(r),func,[]);
end