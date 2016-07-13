%
% PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC converts d'(d-prime) to
% proportion correct for a 2AFC (two-alternative-forced-choice) 
% match-to-sample task, assuming a Differencing Observer model and an 
% unbiased observer
%
% Syntax: [pC]=PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC(dP);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') defined in the range 0<d'<inf 
%
% Example:
%
% [pC]=PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC([0 1 2 3 4 5])
%
% returns:
% 
% pC =
% 
%     0.5000    0.5825    0.7468    0.8765    0.9467    0.9792
%
% The example input argument is an N=6 vector of d' and the output the
% resulting vector of proportion correct
%
% Introduced: Palamedes version 1.0.0 (FK)


function [pC]=PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC(dP)

pC=PAL_ZtoP(dP/sqrt(2)).*PAL_ZtoP(dP/sqrt(6)) + PAL_ZtoP(-1.*dP/sqrt(2)).*PAL_ZtoP(-1.*dP/sqrt(6));