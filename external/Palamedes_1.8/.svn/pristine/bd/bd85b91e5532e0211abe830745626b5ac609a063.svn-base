%PAL_isIdentity     Check whether matrix is an identity matrix
%
%Syntax: eyeIs = PAL_isIdentity(x)
%
%Example 1: PAL_isIdentity([1 0; 0 1]) returns 1
%
%Example 2: PAL_isIdentity([1 0; 1 1]) returns 0
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.4.0 (see History.m)

function eyeIs = PAL_isIdentity(x)

eyeIs = ~isempty(x) && size(x,1) == size(x,2) && sum(sum(x == eye(size(x)))) == size(x,1).*size(x,2);