function localConeToRGCDensityRatioStruct = localConeToRGCDensityRatioFunctionHandle(...
    RGCRFposMicrons, theInputConeMosaic, sampleSizeMicrons)

    
    RGCspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFposMicrons);
    RGCdensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(RGCspacingsMicrons);
    
    coneDensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(theInputConeMosaic.coneRFspacingsMicrons);    
    
    interpolationMethod = 'linear';
    extrapolationMethod = 'linear';
    RGCDensityFunctionHandle = scatteredInterpolant(...
        RGCRFposMicrons(:,1), ...
        RGCRFposMicrons(:,2), ...
        RGCdensities(:), ...
        interpolationMethod, extrapolationMethod);


    coneDensityFunctionHandle = scatteredInterpolant(...
        theInputConeMosaic.coneRFpositionsMicrons(:,1), ...
        theInputConeMosaic.coneRFpositionsMicrons(:,2), ...
        coneDensities(:), ...
        interpolationMethod, extrapolationMethod);

    minX1 = min(RGCRFposMicrons(:,1));
    maxX1 = max(RGCRFposMicrons(:,1));
    minY1 = min(RGCRFposMicrons(:,2));
    maxY1 = max(RGCRFposMicrons(:,2));

    minX2 = min(theInputConeMosaic.coneRFpositionsMicrons(:,1));
    maxX2 = max(theInputConeMosaic.coneRFpositionsMicrons(:,1));
    minY2 = min(theInputConeMosaic.coneRFpositionsMicrons(:,2));
    maxY2 = max(theInputConeMosaic.coneRFpositionsMicrons(:,2));

    minX = max([minX1 minX2]);
    maxX = min([maxX1 maxX2]);
    minY = max([minY1 minY2]);
    maxY = min([maxY1 maxY2]);

    xSamples = ceil((maxX-minX)/sampleSizeMicrons);
    ySamples = ceil((maxY-minY)/sampleSizeMicrons);
    xSupport = linspace(minX, maxX, xSamples);
    ySupport = linspace(minY, maxY, ySamples);

    localConeToRGCDensityRatioStruct = struct(...
        'xSupport', xSupport, ...
        'ySupport', ySupport, ...
        'coneDensityFunctionHandle', coneDensityFunctionHandle, ...
        'RGCDensityFunctionHandle', RGCDensityFunctionHandle ...
        );

end