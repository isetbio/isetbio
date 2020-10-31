function [theConeMosaic, theConeMosaicMetaData] = generateInputConeMosaic(generationMode, eccDegs, sizeDegs,  extraDegsForRGCSurround, coneRFpositionsMicronsInSpatiallyVaryingLattice)
    switch (generationMode)
        case 'equivalent regular hex'
            [theConeMosaic, theConeMosaicMetaData] = generateEquivalentRegularHexConeMosaic(...
                eccDegs, sizeDegs,  extraDegsForRGCSurround, coneRFpositionsMicronsInSpatiallyVaryingLattice);
        otherwise
            error('Unknown generation mode (''%s'') for input cone mosaic', generationMode);
    end
end

function [theConeMosaic, theConeMosaicMetaData] = generateEquivalentRegularHexConeMosaic(eccDegs, sizeDegs,  extraDegsForRGCSurround, coneRFpositionsMicronsInSpatiallyVaryingLattice)
    
    % Determine the mean spacing of the imported cone positions.
    allPairWiseConeDistancesMicrons = pdist2(...
        coneRFpositionsMicronsInSpatiallyVaryingLattice, ...
        coneRFpositionsMicronsInSpatiallyVaryingLattice, ...
    	'euclidean','Smallest',2);
    nearestDistances = squeeze(allPairWiseConeDistancesMicrons(2,:));
    
    % Set the cone spacing in the regular cone mosaic to the mean cone
    % spacing of the imported cone positions
    coneSpacingMicrons = mean(nearestDistances);
    coneApertureMicrons = coneSpacingMicrons * RGCmodels.Watson.constants.coneApertureToDiameterRatio;
    
    % Set the microns-per-degree factor
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));
    micronsPerDegree = sizeMicrons(1)/sizeDegs(1);
    
    % Other mosaic parameters
    coneMosaicResamplingFactor = 5;
    spatialDensity = [0 0.6 0.3 0.1];
    integrationTime = 5/1000;
    sConeFreeRadiusMicrons = 0;
    
    % Generate the regular cone mosaic, with a larger fov than the RGC  mosaic
    % allowing for cones to cover the surrounds of the most distant mRGC   
    theConeMosaic = coneMosaicHex(coneMosaicResamplingFactor, ...
        'fovDegs', sizeDegs + 2 * extraDegsForRGCSurround, ...
        'micronsPerDegree', micronsPerDegree, ...
        'integrationTime', integrationTime, ...
        'customLambda', coneSpacingMicrons, ...
        'customInnerSegmentDiameter', coneApertureMicrons, ...
        'spatialDensity', spatialDensity, ...
        'sConeMinDistanceFactor', 2, ...
        'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );

    % Set the cone quantal absorptions to be appropriate for the sqrt(sum(eccDegs.^2,2))
    fprintf(2,'\nSet cone mosaic quantal efficiencies for the mosaic eccentricity\n');
    fprintf(2,'\nSet macular pigment density for the mosaic eccentricity\n');
    
    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center
    eccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    coneAperturesMicrons = ones(size(conePositionsMicrons,1),1) * coneApertureMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    % Set the metadata struct, which includes the ecc-based cone positions
    % the cone spacings and cone apertures sizes, and the cone types
    % These are all used for connecting cones to the subregions of the 
    % mRGC cell receptive fields.
    theConeMosaicMetaData = struct(...
        'conePositionsMicrons', conePositionsMicrons, ...
        'coneAperturesMicrons', coneAperturesMicrons, ...
        'coneSpacingsMicrons', coneSpacingsMicrons, ...
        'coneTypes', coneTypes ...
        );
end

