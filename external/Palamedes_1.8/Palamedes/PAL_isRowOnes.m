%PAL_isRowOnes     Check whether array is a row vector containing ones
%
%Syntax: rowOnesIs = PAL_isRowOnes(x)
%
%Example 1: PAL_isRowOnes([1 1 1 1]) returns 1
%
%Example 2: PAL_isRowOnes([1 0 1 1]) returns 0
%
%Example 3: PAL_isRowOnes([1 1 1 1]') returns 0
%
%Introduced: Palamedes version 1.0.0 (NP)

function rowOnesIs = PAL_isRowOnes(x)

rowOnesIs = size(x,1) == 1 && sum(ones(size(x,1))==x) == length(x);