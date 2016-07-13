%
% PAL_SDT_2AFC_PCtoDP converts proportion correct to d' (d-prime) for a 
% 2AFC (two-alternative-forced-choice) task
%
% Syntax: [dP]=PAL_SDT_2AFC_DPtoPC(PC);
% 
% returns a scalar, vector or matrix of d' (0<d'<inf) for a scalar, vector 
% or matrix of proportion correct (PC) defined in the range 0<PC<1
%
% Example:
%
% [dP]=PAL_SDT_2AFC_PCtoDP([0.4 .52 0.6 0.85 0.99])
%
% returns:
% 
% dP =
% 
%    -0.3583    0.0709    0.3583    1.4657    3.2900
% 
% The example input consists of a N=5 vector of PC. The output is a
% vector of d's
%
% Introduced: Palamedes version 1.8.0 (FK)

function dP=PAL_SDT_2AFC_PCtoDP(PC)

dP=sqrt(2).*PAL_PtoZ(PC);
