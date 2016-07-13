%
% PAL_SDT_2AFCsameDiff_DPtoPC converts d'(d-prime) to
% proportion correct for a 2AFC (two-alternative-forced-choice) 
% same-different task, assuming a Differencing model and an unbiased 
% observer
%
% Syntax: [pC]=PAL_SDT_2AFCsameDiff_DPtoPC(dP);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') defined in the range 0<d'<inf 
%
% Example:
%
% [pC]=PAL_SDT_2AFCsameDiff_DPtoPC([0 1 2 3 4 5])
% 
% pC =
% 
%     0.5000    0.5733    0.7330    0.8753    0.9555    0.9877
%
% The example input consists of a N=6 vector of d' and the output the
% corresponding N=6 proportion correct
%
% Introduced: Palamedes version 1.0.0 (FK)

function pC=PAL_SDT_2AFCsameDiff_DPtoPC(dP)

pC=PAL_ZtoP(dP./2).^2 + PAL_ZtoP(-1.*dP./2).^2;