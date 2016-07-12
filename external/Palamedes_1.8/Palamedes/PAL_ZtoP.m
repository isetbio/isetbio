%
% PAL_ZtoP converts a z score to a proportion p, where p is the ordinate 
% value of a cumulative normal distribution corresponding to the abscissa 
% value z expressed in standard deviation units
%
% Syntax: [P]=PAL_ZtoP(Z);
%
% returns a scalar, vector or matrix of proportions p ('P') for a scalar,
% vector or matrix of z scores ('Z'), defined in the range -inf<z<inf
%
% Example:
%
% [P]=PAL_ZtoP([-5 0 1.96])
%
% returns:
%
% P =
%
%    0.0000    0.5000    0.9750 
%
% The example input is an N=3 vector of z and the output the resulting
% vector of p
%
% Introduced: Palamedes version 1.0.0 (FK)

function P = PAL_ZtoP(Z)

P = 0.5+0.5*erf((Z)./sqrt(2));