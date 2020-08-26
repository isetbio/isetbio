function [theConeMosaic, theConeMosaicMetaData] = generateConeMosaicForDeconvolution(patchEcc, patchSize, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('sizeUnits', 'degrees', @(x)(ismember(x, {'microns', 'degrees'})));
    p.addParameter('coneDensities', [0.6 0.3 0.1], @(x)((isnumeric(x))&&(numel(x)==3)));
    p.addParameter('mosaicGeometry', 'regular', @(x)(ismember(x, {'regular', 'eccVarying'})));
    p.addParameter('coneMosaicResamplingFactor', 9, @isnumeric);
    p.parse(varargin{:});
    
    sizeUnits = p.Results.sizeUnits;
    mosaicGeometry = p.Results.mosaicGeometry; 
    coneMosaicResamplingFactor = p.Results.coneMosaicResamplingFactor;
    coneDensities = p.Results.coneDensities;
    
    [theConeMosaic, coneMosaicEccDegs, coneMosaicEccMicrons, ...
     coneMosaicSizeDegs, coneMosaicSizeMicrons, extraMicronsForSurroundCones, ...
     conePositionsMicrons, coneSpacingsMicrons, coneAperturesMicrons, coneTypes] = ...
        generateRegularHexMosaicPatch(coneMosaicResamplingFactor, ...
                patchEcc, patchSize, sizeUnits, coneDensities);

    if (strcmp(mosaicGeometry, 'regular'))
        conePositionsDegs = conePositionsMicrons / theConeMosaic.micronsPerDegree;
        coneAperturesDegs = coneAperturesMicrons / theConeMosaic.micronsPerDegree;
    else
        conePositionsDegs = WatsonRGCModel.rhoMMsToDegs(conePositionsMicrons*1e-3);
        coneAperturesDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(coneAperturesMicrons, coneMosaicEccMicrons);
    end
    
    % Assemble metadata struct for the mosaic
    theConeMosaicMetaData = struct(...
        'coneMosaicEccDegs', coneMosaicEccDegs, ...
        'coneMosaicEccMicrons', coneMosaicEccMicrons, ...
        'coneMosaicSizeDegs', coneMosaicSizeDegs, ...
        'coneMosaicSizeMicrons', coneMosaicSizeMicrons, ...
        'conePositionsMicrons', conePositionsMicrons, ...
        'coneAperturesMicrons', coneAperturesMicrons, ...
        'coneSpacingsMicrons', coneSpacingsMicrons, ...
        'conePositionsDegs', conePositionsDegs, ...
        'coneAperturesDegs', coneAperturesDegs, ...
        'coneTypes', coneTypes, ...
        'extraMicronsForSurroundCones', extraMicronsForSurroundCones ...
        );
            
end


function [theConeMosaic, patchEccDegs, patchEccMicrons, ...
      patchSizeDegs, patchSizeMicrons, extraMicronsForSurroundCones, ...
      conePositionsMicrons, coneSpacingsMicrons, coneAperturesMicrons, coneTypes] = generateRegularHexMosaicPatch(...
      coneMosaicResamplingFactor, patchEcc, patchSize, sizeUnits, coneDensities)

    ck = CronerKaplanRGCModel('generateAllFigures', false);
    
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

    patchEccRadius = sqrt(sum(patchEccDegs.^2,2));
    surroundCharacteristicRadiusDegs = ck.surroundRadiusFunction(ck.surroundRadiusParams, patchEccRadius);
    
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
    spatialDensity = [0 coneDensities(1) coneDensities(2) coneDensities(3)];
    sConeFreeRadiusMicrons = 0;
    fprintf('Generating regular hex cone mosaic at eccentricity of (%2.1f,%2.1f) degs,\nwith a FOV of (%2.1f x %2.1f) degs\nand cone densities of %2.2f:%2.2f:%2.2f\n', ...
        patchEccDegs(1), patchEccDegs(2), patchSizeDegs(1), patchSizeDegs(2), coneDensities(1), coneDensities(2), coneDensities(3));
   
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

function [coneSpacingMicrons, coneSpacingMaxMicrons] = medianConeSpacingInPatch(whichEye, patchEccMicrons, patchSizeMicrons)

    w = WatsonRGCModel('generateAllFigures', false);
    posUnits = 'mm'; densityUnits = 'mm^2';
    
    % Cone spacing at the center of the patch
    coneSpacingMicronsAtPatchCenter = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
        patchEccMicrons*1e-3, whichEye, posUnits, densityUnits, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
    % Generate grid of nodes on which to evaluare cone spacing
    xx = 0:coneSpacingMicronsAtPatchCenter:(0.5*patchSizeMicrons(1));
    xx = [-fliplr(xx) xx(2:end)];
    x = patchEccMicrons(1) + xx;
    yy = 0:coneSpacingMicronsAtPatchCenter:(0.5*patchSizeMicrons(2));
    yy = [-fliplr(yy) yy(2:end)];
    y = patchEccMicrons(2) + yy;
    [X,Y] = meshgrid(x,y);
    X = X(:); Y = Y(:);

    % Determine cone spacings at all nodes of the grid
    coneSpacingMicronsList = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
            [X Y]*1e-3, whichEye, posUnits, densityUnits, ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
    % Determine the median cone spacing with the patch
    coneSpacingMicrons = median(coneSpacingMicronsList);
    coneSpacingMaxMicrons = max(coneSpacingMicronsList);
    fprintf('Cone spacing (patch center): %2.1f microns, median over patch: %2.1f\n', ...
        coneSpacingMicronsAtPatchCenter, coneSpacingMicrons);
end
