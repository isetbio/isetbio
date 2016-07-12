%
%PAL_findMax   find value and position of maximum in N-D array
%
%   syntax: [maxim I] = PAL_findMax(array)
%
%   [maxim I] = PAL_findMax(array) returns the value and position of 
%       maximum in an N-D array.
%
%   For vectors, PAL_findMax differs from Matlab's resident 'max' in that 
%   output argument 'I' has two elements: [row, column].
%
%   Example:
%
%       x = zeros(5,5,5);
%       x(1,2,3) = 12;
%       [maxim I] = PAL_findMax(x)
%
%       returns:
%
%       maxim =  
%
%           12
%       I =
%
%           1     2     3
%
% Introduced: Palamedes version 1.1.1 (NP)
% Modified: Palamedes version 1.2.0, 1.6.0, 1.6.3 (see History.m)

function [maxim, Indices] = PAL_findMax(array)

s = size(array);
[maxim, I] = max(reshape(array,[1 prod(s)]));

for dim = length(s):-1:1
    f = prod(s(1:dim-1));     
    Indices(dim) = 1+floor((I-1)/f);
    I = I-(Indices(dim)-1)*f;
end