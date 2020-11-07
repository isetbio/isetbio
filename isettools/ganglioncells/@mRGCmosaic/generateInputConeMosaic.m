function [theConeMosaic, theConeMosaicMetaData] = generateInputConeMosaic(generationMode, eccDegs, sizeDegs,  ...
    extraDegsForRGCSurround, coneRFpositionsMicronsInSpatiallyVaryingLattice, varargin)
    
    % Parse optional params
    p = inputParser();
    p.addParameter('coneMosaicResamplingFactor', 5, @isnumeric);
    p.addParameter('coneMosaicSpatialDensity', [0.6 0.3 0.1], @(x)(isnumeric(x) && numel(x) == 3));
    p.addParameter('coneMosaicIntegrationTime', 5/1000, @isnumeric);
    p.parse(varargin{:});
    coneMosaicResamplingFactor = p.Results.coneMosaicResamplingFactor;
    coneMosaicSpatialDensity = [0 p.Results.coneMosaicSpatialDensity(1) p.Results.coneMosaicSpatialDensity(2) p.Results.coneMosaicSpatialDensity(3)];
    coneMosaicIntegrationTime = p.Results.coneMosaicIntegrationTime;
    switch (generationMode)
        case 'equivalent regular hex'
            [theConeMosaic, theConeMosaicMetaData] = generateEquivalentRegularHexConeMosaic(...
                eccDegs, sizeDegs,  extraDegsForRGCSurround, ...
                coneRFpositionsMicronsInSpatiallyVaryingLattice, ...
                coneMosaicResamplingFactor, ...
                coneMosaicSpatialDensity, ...
                coneMosaicIntegrationTime);
        otherwise
            error('Unknown generation mode (''%s'') for input cone mosaic', generationMode);
    end
end

function [theConeMosaic, theConeMosaicMetaData] = generateEquivalentRegularHexConeMosaic(eccDegs, sizeDegs,  ...
    extraDegsForRGCSurround, coneRFpositionsMicronsInSpatiallyVaryingLattice, ...
    coneMosaicResamplingFactor, coneMosaicSpatialDensity, coneMosaicIntegrationTime)

    % Compute the mosaic ecc radius
    mosaicEccRadiusDegs = sqrt(sum(eccDegs.^2,2));
    
    % Compute cone spacings in the spatially-varying lattice
    eccVaryingConeSpacings = RGCmodels.Watson.convert.positionsToSpacings(coneRFpositionsMicronsInSpatiallyVaryingLattice);
    
    % Compute mean cone spacing within the central region of the spatially-varying lattice
    mosaicCenterPosMicrons = mean(coneRFpositionsMicronsInSpatiallyVaryingLattice,1);
    d = sqrt(sum((bsxfun(@minus,coneRFpositionsMicronsInSpatiallyVaryingLattice, mosaicCenterPosMicrons)).^2,2));
    centralRegionRadius = max(d)/4;
    coneIndicesAtPatchCenter = find(d<=centralRegionRadius);
    
    % Set the cone spacing
    coneSpacingMicronsHexRegMosaic = mean(eccVaryingConeSpacings(coneIndicesAtPatchCenter));
    
    % Set the cone aperture
    coneApertureMicronsHexRegMosaic = coneSpacingMicronsHexRegMosaic * RGCmodels.Watson.constants.coneApertureToDiameterRatio;
    
    % Set the microns-per-degree factor
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, mosaicEccRadiusDegs);
    micronsPerDegree = sizeMicrons(1)/sizeDegs(1);
    
    if (mosaicEccRadiusDegs == 0)
        % Tritanopic area: 0.3 deg in diameter
        tritanopicDiameterDegs = 0.3;
        sConeFreeRadiusMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(0.5*tritanopicDiameterDegs, mosaicEccRadiusDegs);
    else
        sConeFreeRadiusMicrons = 0;
    end
    
    % Generate the regular cone mosaic, with a larger fov than the RGC  mosaic
    % allowing for cones to cover the surrounds of the most distant mRGC   
    theConeMosaic = coneMosaicHex(coneMosaicResamplingFactor, ...
        'fovDegs', sizeDegs + 2 * extraDegsForRGCSurround, ...
        'micronsPerDegree', micronsPerDegree, ...
        'integrationTime', coneMosaicIntegrationTime, ...
        'customLambda', coneSpacingMicronsHexRegMosaic, ...
        'customInnerSegmentDiameter', coneApertureMicronsHexRegMosaic, ...
        'spatialDensity', coneMosaicSpatialDensity, ...
        'sConeMinDistanceFactor', 2, ...
        'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );

    % Set the cone optical density to be appropriate for the outer-segment 
    % length at the mosaic's eccenticity, mosaicEccRadiusDegs
    theConeMosaic.pigment.opticalDensity = eccBasedPhotpigmentOpticalDensities(theConeMosaic, mosaicEccRadiusDegs);
    
    % Set the macular pigment density to be appropriate for the mosaic's
    % eccentricity, mosaicEccRadiusDegs
    theConeMosaic.macular.density = Macular.eccDensity(mosaicEccRadiusDegs);
    
    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center
    eccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccMicrons);
    conePositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(conePositionsMicrons*1e-3);
    
    % Find indices of cones that are not to be connected to the RGC centers
    % These are cones positioned outside the FOV of the RGC mosaic
    d = bsxfun(@minus,abs(cmStruct.coneLocsMicrons), sizeMicrons/2);
    idxConesOutsideRGCmosaic = find( (d(:,1) > 0) | (d(:,2) > 0));
    
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    % Set the metadata struct, which includes the ecc-based cone positions
    % and the cone types. These are all used for connecting cones to the subregions of the 
    % mRGC cell receptive fields.
    theConeMosaicMetaData = struct(...
        'indicesOfConesNotConnectingToRGCcenters', idxConesOutsideRGCmosaic, ...
        'conePositionsMicrons', conePositionsMicrons, ...
        'conePositionsDegs', conePositionsDegs, ...
        'coneSpacingsMicrons', RGCmodels.Watson.convert.positionsToSpacings(conePositionsMicrons), ...
        'coneSpacingsDegs', RGCmodels.Watson.convert.positionsToSpacings(conePositionsDegs), ...
        'coneTypes', coneTypes ...
        );
end

function peripheralOpticalDensities = eccBasedPhotpigmentOpticalDensities(theConeMosaic, mosaicEccDegs) 
    % Scanned data (eccentricity, osLength) from Figure 1 (right panel)
    % Banks, Sekuler and Anderson (1991). Peripher spatial vision: limits
    % imposed by optics, photoreceptors and receptor pooling
    eccDegsBanksSekulerAnderson = [0 1.82 4.86 9.86 19.78 39.9];
    osLengthMicronsBanksSekulerAnderson = [47.81 26.16 21.2 21.2 21.2 13.22];
    
    % Optical density = specific density x outer segment length (microns)
    % Compute specific densities assuming cone foveal length that of Banks,
    % Sekuler and Andrerson
    fovealOpticalDensities = theConeMosaic.pigment.opticalDensity;
    fovealOuterSegmentLengthMicrons = osLengthMicronsBanksSekulerAnderson(1);
    specificDensities = fovealOpticalDensities / fovealOuterSegmentLengthMicrons;
    
    % Given these specific densities, compute peripheral optical density
    % at the desired eccentricity;
    peripheralOuterSegmentLengthMicrons = interp1(eccDegsBanksSekulerAnderson, ...
        osLengthMicronsBanksSekulerAnderson, mosaicEccDegs, 'phip');
    peripheralOpticalDensities = specificDensities * peripheralOuterSegmentLengthMicrons;
end