%
% PAL_SDT_1AFCsameDiff_IndMod_DPtoPHF converts d' (d-prime) with 
% criterion C into proportion hits and proportion false alarms for a 1AFC 
% (one-alternative-forced-choice) same-different task under the 
% Independent Observation model
%
% Syntax: [pHF]=PAL_SDT_1AFCsameDiff_IndMod_DPtoPHF(dP,C);
%
% returns a Nx2 matrix of N proportion hits and proportion false alarms 
% ('pHF') for a scalar or N-length vector of d' ('dP') and criterion C 
% ('C'), defined in the ranges 0<d'<inf and -inf<C<inf
%
% Example:
%
% [pHF]=PAL_SDT_1AFCsameDiff_IndMod_DPtoPHF([.5 .6 .7],[-1.0 0 1.0])
%
% returns:
%
% pHF =
%
%    0.8529    0.8292
%    0.5278    0.4722
%    0.1825    0.1370
%
% The input arguments are two N=3 vectors of d' and C. The first column 
% in the 3 x 2 matrix output gives the proportion of hits and the second      
% column the proportion of false alarms
%
% Introduced: Palamedes version 1.0.0 (FK)

function [pHF]=PAL_SDT_1AFCsameDiff_IndMod_DPtoPHF(dP,C)

PCmax=PAL_ZtoP(dP/2).^2 + PAL_ZtoP(-1.*dP/2).^2;
zF=-C-PAL_PtoZ(PCmax);
zH=-zF-2.*C;
pHF(:,1)=PAL_ZtoP(zH);
pHF(:,2)=PAL_ZtoP(zF);