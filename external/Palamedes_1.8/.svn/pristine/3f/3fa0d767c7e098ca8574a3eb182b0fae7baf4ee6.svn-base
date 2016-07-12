%
%PAL_isOrthogonal  Determines whether a set of vectors is orthogonal
%
%   syntax: [ortho nonOrtho] = PAL_isOrthogonal(M)
%
%   ortho: logical 1 if row vectors are all mutually orthogonal, 0 if not.
%   nonOrtho: array listing nor orthogonal pairings of row vectors.
%
%Internal function
%
%Introduced: Palamedes version 1.6.2 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [ortho, nonOrtho] = PAL_isOrthogonal(M)

ortho = logical(1);
nonOrtho = [];
count = 0;

for C1 = 1:size(M,1)-1
    for C2 = C1+1:size(M,1)
        dotProd = sum(M(C1,:).*M(C2,:));
        if dotProd            
            ortho = 0;
            count = count + 1;
            nonOrtho(count,:) = [C1 C2 dotProd];
        end
    end
end