%
% PAL_SDT_2AFCmatchSample_IndMod_DPtoPC converts d'(d-prime) to
% proportion correct for a 2AFC (two-alternative-forced-choice) 
% match-to-sample task, assuming an Independent Observer model and an 
% unbiased observer
%
% Syntax: [pC]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPC(dP);
% 
% returns a matrix of proportion correct ('pC') for a matrix of d' ('dP') 
% defined in the range 0<d'<inf 
%
% Example:
%
% [pC]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPC([0 1 2 3 4 5])
%
% returns:
% 
% pC =
% 
%     0.5000    0.5997    0.7877    0.9185    0.9750    0.9936
%
% 
% The example input consists of a N=6 vector of d' and the output the
% corresponding N=6 proportion correct
%
% Introduced: Palamedes version 1.0.0 (FK)

function [pC]=PAL_SDT_2AFCmatchSample_IndMod_DPtoPC(dP)

pC=PAL_ZtoP(dP/sqrt(2)).*PAL_ZtoP(dP/2) + PAL_ZtoP(-1.*dP/sqrt(2)).*PAL_ZtoP(-1.*dP/2);