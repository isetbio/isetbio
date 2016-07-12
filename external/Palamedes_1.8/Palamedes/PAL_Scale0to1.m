%PAL_Scale0to1      Rescale N-D array linearly such that highest value in 
%   array equals 1 and lowest value in array equals 0.
%
%   syntax: y = PAL_Scale0to1(x)
%
%Example:
%
%   y = PAL_Scale0to1([1 2 3; 3 4 5])
%
%   y = 
%       0       0.2500  0.5000
%       0.5000  0.7500  1.0000
%
%Introduced: Palamedes version 1.0.0 (FK)
%Modified: Palamedes version 1.2.0 (see History.m)

function x = PAL_Scale0to1(x)

minx = min(min(x));
maxx = max(max(x));

while ~isscalar(minx)
    minx = min(minx);
    maxx = max(maxx);
end

x=(x-minx)./(maxx-minx);
    