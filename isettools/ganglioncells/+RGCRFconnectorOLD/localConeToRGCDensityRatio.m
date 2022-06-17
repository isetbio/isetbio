function localConeToRGCDensityRatio = localConeToRGCDensityRatio(RGCRFposMicrons, theInputConeMosaic)

    % Compute the local RGCRF spacings for all RGCs in the mosaic from
    % their positions
    RGCspacings = RGCmodels.Watson.convert.positionsToSpacings(RGCRFposMicrons);

    % Compute the corresponding local RGC density
    RGCdensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(RGCspacings);

    % Compute the cone densities for all cone positions in the mosaic
    coneSpacings = theInputConeMosaic.coneRFspacingsMicrons;
    coneDensities = RGCmodels.Watson.convert.spacingToDensityForHexGrid(coneSpacings);

    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;

    % Find the closest RGC to each cone
    [~,closestConeIndices] = pdist2(allConePositions, RGCRFposMicrons, '', 'smallest', 1);

    % Preallocate memory
    rgcsNum = size(RGCRFposMicrons,1);
    localConeToRGCDensityRatio = zeros(1, rgcsNum);

    % Do it
    parfor iRGC = 1:rgcsNum
        theNearestConeIndex = closestConeIndices(iRGC);
        localConeToRGCDensityRatio(iRGC) = coneDensities(theNearestConeIndex) / RGCdensities(iRGC);
    end
end