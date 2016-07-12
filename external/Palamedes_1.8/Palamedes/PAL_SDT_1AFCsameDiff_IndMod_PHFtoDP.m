%
% PAL_SDT_1AFCsameDiff_IndMod_PHFtoDP converts proportion hits and 
% proportion false alarms into d'(d-prime) and criterion k for a 1AFC 
% (one-alternative-forced-choice) same-different task under the
% Independent Observation model
%
% Syntax: [dP C]=PAL_SDT_1AFCsameDiff_IndMod_PHFtoDP(pHF);
% 
% returns a scalar or N-length vector of d' ('dP') and criterion C ('C') 
% for an Nx2 input matrix of N proportion hits and proportion false alarms 
% ('pHF') defined in the range 0<p<1 (p=proportion)
%
% Example:
%
% [dP C]=PAL_SDT_1AFCsameDiff_IndMod_PHFtoDP([0.6 0.1; 0.8 0.1; 0.9 0.6])
%
% returns:
% 
% dP =
% 
%     2.2835
%     2.8342
%     1.7808
% 
% 
% C =
% 
%     0.5141
%     0.2200
%    -0.7674
%
% The example input argument is a 3 x 2 matrix, with each row (demarcated
% by a semi-colon) consisting of a proportion of hits and proportion of 
% false alarms.  The columns in the output are the resulting N=3 
% vectors of dP and C
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.6.3 (see History.m)

function [dP, C]=PAL_SDT_1AFCsameDiff_IndMod_PHFtoDP(pHF)

zH=PAL_PtoZ(pHF(:,1));
zF=PAL_PtoZ(pHF(:,2));

diff=(zH-zF)./2;
PCmax=PAL_ZtoP(diff);

val=0.5*(1+sqrt(2.*PCmax-1));
dP=2*PAL_PtoZ(val);
C=-0.5.*(zH+zF);