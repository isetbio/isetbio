function [spacings, idx] = positionsToSpacings(positions)
% positionsToSpacings computes the median of the 6 smallest distances
% between points in a given Nx2 matrix of positions.
%
% (Copilot comments).
%
% Inputs:
%   positions - An Nx2 matrix where each row represents the (x, y)
%               coordinates of a point.
%
% Outputs:
%   spacings - A 1x2 vector containing the median of the smallest
%              distances for each dimension.
%   idx - Indices of the 6 nearest neighbors for each point.
%
% Example:
%   positions = [1, 2; 3, 4; 5, 6];
%   [spacings, idx] = positionsToSpacings(positions);

assert(size(positions,2) == 2, 'Positions must be an Nx2 matrix');
[allPairWiseDistances, idx] = pdist2(positions, positions, ...
    'euclidean','Smallest',7);

% Mean spacing removing the first point which is 0 (distance to itself)
spacings = median(allPairWiseDistances(2:end,:),1);

end
