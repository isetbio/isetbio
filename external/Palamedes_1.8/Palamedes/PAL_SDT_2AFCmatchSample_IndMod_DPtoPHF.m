%
% PAL_SDT_2AFCmatchSample_IndMod_DPtoPHF converts d' (d-prime) with 
% criterion C into proportion hits and proportion false alarms for a 2AFC 
% (two-alternative-forced-choice) match-to-sample task under the 
% Differencing model
%
% Syntax: [pHF]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPHF(dP,C);
%
% returns a Nx2 matrix of N proportion hits and proportion false alarms 
% ('pHF') for a scalar or N-length vector of d' ('dP') and criterion C 
% ('C'), defined in the ranges 0<d'<inf and -inf<C<inf
%
% Example:
%
% [pHF]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPHF([0 2.5 3],[-2.5 0 2.5])
%
% returns:
% 
% pHF =
% 
%     0.9938    0.9938
%     0.8639    0.1361
%     0.1346    0.0000
%
% The input arguments are two N=3 vectors of d' and C. The first column 
% in the 3 x 2 matrix output gives the proportion of hits and the second      
% column the proportion of false alarms
%
% Introduced: Palamedes version 1.0.0 (FK)

function [pHF]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPHF(dP,C)

PCmax=PAL_SDT_2AFCmatchSample_IndMod_DPtoPC(dP);

zF=-C-PAL_PtoZ(PCmax);
zH=-zF-2.*C;
pHF(:,1)=PAL_ZtoP(zH);
pHF(:,2)=PAL_ZtoP(zF);