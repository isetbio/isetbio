function densityRatiosMap = coneToRGCDensityRatiosMap(obj)
% Employ the obj.coneToRGCDensityRatioComputeStruct to compute the
% local cone-to-RGC density ratios at the current RGC RF positions

    % Compute the spatial support for the map
    X = obj.RGCRFpositionsMicrons(:,1);
    Y = obj.RGCRFpositionsMicrons(:,2);

    % Compute the density map for the cones
    coneDensityMap = obj.coneToRGCDensityRatioComputeStruct.coneDensityFunctionHandle(X,Y);

    % Compute the density map for the RGCs
    rgcDensityMap = obj.coneToRGCDensityRatioComputeStruct.RGCDensityFunctionHandle(X,Y);

    % Compute the ratio of densities map
    densityRatiosMap = coneDensityMap ./ rgcDensityMap;

end