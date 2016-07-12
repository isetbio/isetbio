%
% PAL_PtoZ converts a proportion p to a z score, where z is the value on 
% the abscissa of a cumulative normal distribution in units of standard 
% deviation and p the corresponding ordinate value
%
% Syntax: [Z]=PAL_PtoZ(P);
%
% returns a scalar, vector or matrix of z ('Z') for a scalar, 
% vector or matrix of proportion p ('P'), defined in the range 0<p<1
%
% Example:
%
% [Z]=PAL_PtoZ([0.0 0.5 0.975])
%
% returns:
%
% Z =
%
%      -Inf         0    1.9600
%
% The example input is an N=3 vector of p and the resulting output a 
% vector of z
%
% Introduced: Palamedes version 1.0.0 (FK)

function Z = PAL_PtoZ(P)

Z = sqrt(2)*erfinv(P.*2-1);