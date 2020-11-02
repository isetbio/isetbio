function spacings = positionsToSpacings(positions)
   % Compute the smallest spacing for each RF from the positions of all RFs in a mosaic
   assert(size(positions,2) == 2, 'Positions must be an Nx2 matrix');
   allPairWiseDistances = pdist2(positions, positions, ...
    	'euclidean','Smallest',2);
   spacings = squeeze(allPairWiseDistances(2,:));
end

