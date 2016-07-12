%
%PAL_nansum  Emulates Matlab's nansum
%
%   Included in palamedes for compatibility with core GNU Octave
%
%   syntax: x = PAL_nansum(x)
%
%   Internal function
%
%Introduced: Palamedes version 1.4.1 (NP)

function [ x ] = PAL_nansum( x )

x = sum(x(~isnan(x)));

end

