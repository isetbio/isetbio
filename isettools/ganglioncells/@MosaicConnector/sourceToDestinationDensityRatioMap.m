function densityRatioMap = sourceToDestinationDensityRatioMap(obj)
% Employ the obj.sourceToDestinationDensityRatioComputeStruct to compute the
% local source:destination density ratios at the current positions of the
% destination lattice RF positions

    % Compute the spatial support for the map
    X = obj.destinationLattice.RFpositionsMicrons(:,1);
    Y = obj.destinationLattice.RFpositionsMicrons(:,2);

    % Compute the source lattice density map
    sourceDensityMap = obj.sourceToDestinationDensityRatioComputeStruct.sourceDensityFunctionHandle(X,Y);

    % Compute the destination lattice density map
    destinationDensityMap = obj.sourceToDestinationDensityRatioComputeStruct.destinationDensityFunctionHandle(X,Y);

    % Compute the ratio of densities map
    densityRatioMap = sourceDensityMap ./ destinationDensityMap;

end