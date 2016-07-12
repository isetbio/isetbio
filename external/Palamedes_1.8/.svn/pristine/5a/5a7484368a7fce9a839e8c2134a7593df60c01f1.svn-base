%
% PAL_SDT_2AFCsameDiff_PCtoDP converts proportion correct to 
% d'(d-prime) for a 2AFC (two-alternative-forced-choice) same-different
% task, assuming a Differencing Observer model and an unbiased 
% observer
%
% Syntax: [dP]=PAL_SDT_2AFCsameDiff_PCtoDP(pC);
% 
% returns a scalar, vector or matrix of d' ('dP') for an input scalar, 
% vector or matrix of proportion correct ('pC'), defined in the range 
% 0<p<1 (p=proportion)
%
% Example:
%
% dP=PAL_SDT_2AFCsameDiff_PCtoDP([.5 .6 .7 .8 .9])
%
% returns:
% 
% dP =
% 
%          0    1.1872    1.8022    2.4246    3.2368
%
% The example input is an N=5 vector of proportion correct and the output 
% a vector of d's
%
% Introduced: Palamedes version 1.0.0 (FK)

function dP=PAL_SDT_2AFCsameDiff_PCtoDP(pC)

val=0.5.*(1+(2.*pC-1).^0.5);
dP=2*PAL_PtoZ(val);