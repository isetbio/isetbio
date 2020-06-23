function runPhaseX(runParams)
    
    % Compute cone spacing in microns for the retinal position corresponding 
    % to the center of the mosaic and the right eye
    coneMosaicCenterPositionMM = runParams.rgcMosaicPatchPosMicrons * 1e-3;
    whichEye = 'right';
    
    w = WatsonRGCModel('generateAllFigures', false);
    posUnits = 'mm'; densityUnits = 'mm^2';
    coneSpacingMicrons = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
        coneMosaicCenterPositionMM, whichEye, posUnits, densityUnits, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
    % Compute the cone mosaic FOV in degrees
    extraMicronsForSurroundCones = 200;
    regHexConeMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons + 2*extraMicronsForSurroundCones*[1 1];
    coneMosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMM);
    fovDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(regHexConeMosaicPatchSizeMicrons, coneMosaicEccDegs);
    
    resamplingFactor = 5;
    spatialDensity = [0 0.5 0.25 0.15];
    sConeFreeRadiusMicrons = 0;
    
    theConeMosaic = coneMosaicHex(resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'integrationTime', 5/1000, ...
        'customLambda', coneSpacingMicrons, ...
        'customInnerSegmentDiameter', coneSpacingMicrons * 0.7, ...
        'spatialDensity', spatialDensity, ...
        'sConeMinDistanceFactor', 2, ...
        'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );

    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, runParams.rgcMosaicPatchPosMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    
    % Load mRGC RF data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');

  	mRGCRFroi.center = runParams.rgcMosaicPatchPosMicrons;
    mRGCRFroi.size = regHexConeMosaicPatchSizeMicrons - 2.0*extraMicronsForSurroundCones;
    
    % Crop midget mosaic
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
        cropRGCmosaic(RGCRFPositionsMicrons, RGCRFSpacingsMicrons,  desiredConesToRGCratios, mRGCRFroi);
    
    extraMicrons = 20;
    xyRange = max([
        max(conePositionsMicrons(:,1))-min(conePositionsMicrons(:,1))
        max(conePositionsMicrons(:,2))-min(conePositionsMicrons(:,2))]);
    xyRange = 0.5*xyRange + extraMicrons;
    
    
    xRange = mean(RGCRFPositionsMicrons(:,1)) + xyRange*[-1 1];
    yRange = mean(RGCRFPositionsMicrons(:,2)) + xyRange*[-1 1];
    
    xx = cosd(0:10:360);
    yy = sind(0:10:360);
    figure(1); clf;
    subplot(1,2,1);
    hold on
    for k = 1:size(conePositionsMicrons,1)
        r = 0.5*coneSpacingsMicrons(k);
        plot(conePositionsMicrons(k,1)+r*xx, conePositionsMicrons(k,2)+r*yy, 'r-');
    end
    axis 'equal'; axis 'square'
    set(gca, 'XLim', xRange, 'YLim', yRange);
   
    
    subplot(1,2,2);
    hold on;
    for k = 1:size(RGCRFPositionsMicrons,1)
        r = 0.5*RGCRFSpacingsMicrons(k);
        plot(RGCRFPositionsMicrons(k,1)+r*xx, RGCRFPositionsMicrons(k,2)+r*yy, 'k-');
    end
    axis 'equal'; axis 'square'
    set(gca, 'XLim', xRange, 'YLim', yRange);
   
end


function [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
    cropRGCmosaic(RGCRFPositionsMicrons, RGCRFSpacingsMicrons,  desiredConesToRGCratios, roi)

    % Find RGCs within the roi
    idxRGC = positionsWithinROI(roi, RGCRFPositionsMicrons);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
end

function indices = positionsWithinROI(roi, positions)
    d = bsxfun(@minus,positions, roi.center);
    ecc = sqrt(sum(positions.^2,2));
    indices = find((abs(d(:,1)) <= 0.5*roi.size(1)) & (abs(d(:,2)) <= 0.5*roi.size(2)));
    if (isempty(indices))
        d = sqrt(sum(d.^2,2));
        [~,indices] = min(d);
    else
        % Re-order according to increasing eccentricity
        [~,sortedIdx] = sort(ecc(indices), 'ascend');
        indices = indices(sortedIdx);
    end
end
