function [theConeMosaic, theConeMosaicMetaData] = generateConeMosaicForDeconvolution(patchEcc, patchSize, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('sizeUnits', 'degrees', @(x)(ismember(x, {'microns', 'degrees'})));
    p.addParameter('mosaicGeometry', 'regular', @(x)(ismember(x, {'regular', 'eccVarying'})));
    p.addParameter('coneMosaicResamplingFactor', 9, @isnumeric);
    p.parse(varargin{:});
    
    sizeUnits = p.Results.sizeUnits;
    mosaicGeometry = p.Results.mosaicGeometry; 
    coneMosaicResamplingFactor = p.Results.coneMosaicResamplingFactor;
    
    [theConeMosaic, coneMosaicEccDegs, coneMosaicEccMicrons, ...
     coneMosaicSizeDegs, coneMosaicSizeMicrons, extraMicronsForSurroundCones, ...
     conePositionsMicrons, coneSpacingsMicrons, coneAperturesMicrons, coneTypes] = generateRegularHexMosaicPatch(coneMosaicResamplingFactor, patchEcc, patchSize, sizeUnits);

    % Assemble metadata struct for the mosaic
    theConeMosaicMetaData = struct(...
        'coneMosaicEccDegs', coneMosaicEccDegs, ...
        'coneMosaicEccMicrons', coneMosaicEccMicrons, ...
        'coneMosaicSizeDegs', coneMosaicSizeDegs, ...
        'coneMosaicSizeMicrons', coneMosaicSizeMicrons, ...
        'conePositionsMicrons', conePositionsMicrons, ...
        'coneSpacingsMicrons', coneSpacingsMicrons, ...
        'coneAperturesMicrons', coneAperturesMicrons, ...
        'coneTypes', coneTypes, ...
        'extraMicronsForSurroundCones', extraMicronsForSurroundCones ...
        );
            
end


function [theConeMosaic, patchEccDegs, patchEccMicrons, ...
      patchSizeDegs, patchSizeMicrons, extraMicronsForSurroundCones, ...
      conePositionsMicrons, coneSpacingsMicrons, coneAperturesMicrons, coneTypes] = generateRegularHexMosaicPatch(...
      coneMosaicResamplingFactor, patchEcc, patchSize, sizeUnits)

    ck = CronerKaplanRGCModel('generateAllFigures', true);
    
    if (strcmp(sizeUnits, 'degrees'))
        % patch eccentricity in degs and microns
        patchEccDegs = patchEcc;
        patchEccMicrons = 1e3 * WatsonRGCModel.rhoDegsToMMs(patchEccDegs);
        
        % patch size in degs
        patchSizeDegs = patchSize;
    else
        % patch eccentricity in degs and microns
        patchEccMicrons = patchEcc;
        patchEccDegs = WatsonRGCModel.rhoMMsToDegs(patchEccMicrons*1e-3);
        
        % patch size in degs
        patchSizeDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(patchSize, sqrt(sum(patchEccMicrons.^2,2)));
    end

    surroundCharacteristicRadiusDegs = ck.surroundRadiusFunction(ck.surroundRadiusParams, sqrt(sum(patchEccDegs.^2,2)));
    
    % add extra degs to each side of the mosaic
    extraDegsForSurroundCones = 2*3*surroundCharacteristicRadiusDegs;
    patchSizeDegs = patchSizeDegs+extraDegsForSurroundCones;
    % patch size in microns
    patchSizeMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(patchSizeDegs, sqrt(sum(patchEccDegs.^2,2)));
    extraMicronsForSurroundCones = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(extraDegsForSurroundCones, sqrt(sum(patchEccDegs.^2,2)));
    
    % Compute microns per degree on the retina
    micronsPerDegree = patchSizeMicrons(2)/patchSizeDegs(2);
    
    % Determine the median cone spacing with the patch
    whichEye = 'right';
    coneSpacingMicrons = medianConeSpacingInPatch(whichEye, patchEccMicrons, patchSizeMicrons);
    coneApertureMicrons = coneSpacingMicrons * WatsonRGCModel.coneApertureToDiameterRatio;
    
    % Generate reg hex cone mosaic
    spatialDensity = [0 0.6 0.3 0.1];
    sConeFreeRadiusMicrons = 0;
    fprintf('Generating regular hex cone mosaic with FOV %2.1f x %2.1f degs\n', patchEccDegs(1), patchEccDegs(2));
    
    theConeMosaic = coneMosaicHex(coneMosaicResamplingFactor, ...
        'fovDegs', patchSizeDegs, ...
        'micronsPerDegree', micronsPerDegree, ...
        'integrationTime', 5/1000, ...
        'customLambda', coneSpacingMicrons, ...
        'customInnerSegmentDiameter', coneApertureMicrons, ...
        'spatialDensity', spatialDensity, ...
        'sConeMinDistanceFactor', 2, ...
        'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );

    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, patchEccMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    coneAperturesMicrons = ones(size(conePositionsMicrons,1),1) * coneApertureMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
end