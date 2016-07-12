%
% PAL_SDT_2AFC_DPtoPC converts d'(d-prime) into proportion correct for a 
% 2AFC (two-alternative-forced-choice) task
%
% Syntax: [PC]=PAL_SDT_2AFC_DPtoPC(dP);
% 
% returns a scalar, vector or matrix of proportion correct for a scalar, 
% vector or matrix of d' ('dP') defined in the ranges 0<d'<inf
%
% Example:
%
% [PC]=PAL_SDT_2AFC_DPtoPC([.5 1 1.5 2 2.5])
%
% returns:
% 
% PC =
% 
%     0.6382    0.7602    0.8556    0.9214    0.9615
% 
% The example input consists of a N=5 vector of d'. The output is a
% vector of PCs
%
% Introduced: Palamedes version 1.8.0 (FK)

function PC=PAL_SDT_2AFC_DPtoPC(dP)

PC=PAL_ZtoP(dP./sqrt(2));