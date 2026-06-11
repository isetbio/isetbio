function [spacings, idx] = positionsToSpacings(positions)
   % Compute the 6 smallest distances of each RF in a mosaic
   assert(size(positions,2) == 2, 'Positions must be an Nx2 matrix');
   [allPairWiseDistances, idx] = pdist2(positions, positions, ...
    	'euclidean','Smallest',7);
   % Mean spacing removing the first point which is 0 (distance to itself)
   spacings = median(allPairWiseDistances(2:end,:),1);
end

