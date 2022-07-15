function computeConeToRGCDensityRatioComputeStruct(obj, RGCRFposMicrons, RGCRFspacingsMicrons, samplingIntervalMicrons)
 
    
    if (isempty(RGCRFspacingsMicrons))
        % Compute RGCRF spacings from their initial positions
        obj.RGCRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFposMicrons);
    else
        % Use supplied spacings
        assert(numel(RGCRFspacingsMicrons) == size(RGCRFposMicrons,1), ...
            'size of RGCRFspacings (%d) is incosistent with RGCRFposMicrons (%d)', numel(RGCRFspacingsMicrons), size(RGCRFposMicrons,1));
        
        obj.RGCRFspacingsMicrons = reshape(RGCRFspacingsMicrons, [1 size(RGCRFposMicrons,1)]);
    end

    % Compute RGC densities from their spacings
    RGCdensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(obj.RGCRFspacingsMicrons);
    
    % Compute cone densities from their positions
    coneDensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(obj.inputConeMosaic.coneRFspacingsMicrons);    
    
    interpolationMethod = 'linear';
    extrapolationMethod = 'linear';
    RGCDensityFunctionHandle = scatteredInterpolant(...
        RGCRFposMicrons(:,1), ...
        RGCRFposMicrons(:,2), ...
        RGCdensities(:), ...
        interpolationMethod, extrapolationMethod);


    coneDensityFunctionHandle = scatteredInterpolant(...
        obj.inputConeMosaic.coneRFpositionsMicrons(:,1), ...
        obj.inputConeMosaic.coneRFpositionsMicrons(:,2), ...
        coneDensities(:), ...
        interpolationMethod, extrapolationMethod);

    minX1 = min(RGCRFposMicrons(:,1));
    maxX1 = max(RGCRFposMicrons(:,1));
    minY1 = min(RGCRFposMicrons(:,2));
    maxY1 = max(RGCRFposMicrons(:,2));

    minX2 = min(obj.inputConeMosaic.coneRFpositionsMicrons(:,1));
    maxX2 = max(obj.inputConeMosaic.coneRFpositionsMicrons(:,1));
    minY2 = min(obj.inputConeMosaic.coneRFpositionsMicrons(:,2));
    maxY2 = max(obj.inputConeMosaic.coneRFpositionsMicrons(:,2));

    minX = max([minX1 minX2]);
    maxX = min([maxX1 maxX2]);
    minY = max([minY1 minY2]);
    maxY = min([maxY1 maxY2]);

    xSamples = ceil((maxX-minX)/samplingIntervalMicrons);
    ySamples = ceil((maxY-minY)/samplingIntervalMicrons);
    xSupport = linspace(minX, maxX, xSamples);
    ySupport = linspace(minY, maxY, ySamples);

    obj.coneToRGCDensityRatioComputeStruct = struct(...
        'xSupport', xSupport, ...
        'ySupport', ySupport, ...
        'coneDensityFunctionHandle', coneDensityFunctionHandle, ...
        'RGCDensityFunctionHandle', RGCDensityFunctionHandle ...
        );

end