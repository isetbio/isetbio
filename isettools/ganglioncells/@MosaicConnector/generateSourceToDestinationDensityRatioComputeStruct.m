function generateSourceToDestinationDensityRatioComputeStruct(obj, samplingIntervalMicrons)

    % Compute source lattice density from the spacings
    sourceDensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(obj.sourceLattice.RFspacingsMicrons);
    
    % Compute destination lattice density from the spacings
    destinationDensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(obj.destinationLattice.RFspacingsMicrons);
    
    
    interpolationMethod = 'linear';
    extrapolationMethod = 'linear';

    destinationLatticeDensityFunctionHandle = scatteredInterpolant(...
        obj.destinationLattice.RFpositionsMicrons(:,1), ...
        obj.destinationLattice.RFpositionsMicrons(:,2), ...
        destinationDensities(:), ...
        interpolationMethod, extrapolationMethod);

    sourceLatticeDensityFunctionHandle = scatteredInterpolant(...
        obj.sourceLattice.RFpositionsMicrons(:,1), ...
        obj.sourceLattice.RFpositionsMicrons(:,2), ...
        sourceDensities(:), ...
        interpolationMethod, extrapolationMethod);

    minX1 = min(obj.destinationLattice.RFpositionsMicrons(:,1));
    maxX1 = max(obj.destinationLattice.RFpositionsMicrons(:,1));
    minY1 = min(obj.destinationLattice.RFpositionsMicrons(:,2));
    maxY1 = max(obj.destinationLattice.RFpositionsMicrons(:,2));

    minX2 = min(obj.sourceLattice.RFpositionsMicrons(:,1));
    maxX2 = max(obj.sourceLattice.RFpositionsMicrons(:,1));
    minY2 = min(obj.sourceLattice.RFpositionsMicrons(:,2));
    maxY2 = max(obj.sourceLattice.RFpositionsMicrons(:,2));

    minX = max([minX1 minX2]);
    maxX = min([maxX1 maxX2]);
    minY = max([minY1 minY2]);
    maxY = min([maxY1 maxY2]);

    xSamples = ceil((maxX-minX)/samplingIntervalMicrons);
    ySamples = ceil((maxY-minY)/samplingIntervalMicrons);

    obj.sourceToDestinationDensityRatioComputeStruct = struct(...
        'xSupport', linspace(minX, maxX, xSamples), ...
        'ySupport', linspace(minY, maxY, ySamples), ...
        'sourceDensityFunctionHandle', sourceLatticeDensityFunctionHandle, ...
        'destinationDensityFunctionHandle', destinationLatticeDensityFunctionHandle ...
        );

end
