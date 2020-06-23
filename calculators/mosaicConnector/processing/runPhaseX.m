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
    coneMosaicEccDegs = sqrt(sum((WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMM)).^2,2))
    fovDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(regHexConeMosaicPatchSizeMicrons, coneMosaicEccDegs)
    
    
    % Generate reg hex cone mosaic
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
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, runParams.rgcMosaicPatchPosMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    % Load mRGC RF data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');

    % Crop midget mosaic to the size and position of the cone mosaic, leaving enough space for the surround cones
  	mRGCRFroi.center = runParams.rgcMosaicPatchPosMicrons;
    mRGCRFroi.size = regHexConeMosaicPatchSizeMicrons - 2.0*extraMicronsForSurroundCones;
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
        cropRGCmosaic(RGCRFPositionsMicrons, RGCRFSpacingsMicrons,  desiredConesToRGCratios, mRGCRFroi);
    
    % Visualize mosaics to be connected
    visualizeMosaics(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, extraMicronsForSurroundCones, coneMosaicCenterPositionMM)
end

function visualizeMosaics(conePositionsMicrons, coneSpacingsMicrons, coneTypes, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, extraMicronsForSurroundCones, coneMosaicCenterPositionMM)
    
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID 
    
    extraMicrons = 20;
    xyRange = max([
        max(conePositionsMicrons(:,1))-min(conePositionsMicrons(:,1))
        max(conePositionsMicrons(:,2))-min(conePositionsMicrons(:,2))]);
    xyRange = 0.5*xyRange + extraMicrons;
    
    
    xRange = mean(RGCRFPositionsMicrons(:,1)) + xyRange*[-1 1];
    yRange = mean(RGCRFPositionsMicrons(:,2)) + xyRange*[-1 1];
    
   
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'figureWidthInches', 30, ...
            'figureHeightInches', 15);
        
    hFig = figure(1); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 1, ...
        'colsNum', 2, ...
        'leftMargin', 0.04, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.07, ...
        'bottomMargin', 0.07, ...
        'rightMargin', 0.00, ...
        'topMargin', 0.01);
    
    
    xx = cosd(0:10:360);
    yy = sind(0:10:360);
    
    
    ax = theAxesGrid{1,1};
    hold(ax, 'on');
    for k = 1:size(conePositionsMicrons,1)
        r = 0.5*coneSpacingsMicrons(k);
        switch (coneTypes(k))
            case LCONE_ID 
                colorRGB = [1 0 0];
            case MCONE_ID
                colorRGB = [0 0.6 0];
            case SCONE_ID
                colorRGB = [0 0 1];
        end
        patch(ax, conePositionsMicrons(k,1)+r*xx, conePositionsMicrons(k,2)+r*yy,  colorRGB);
    end
    axis(ax, 'equal'); axis(ax, 'square')
    set(ax, 'XLim', xRange, 'YLim', yRange);
    xlabel(ax, 'microns');
    mosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMM);
    title(ax,sprintf('reg hex cone mosaic patch (x,y) = (%2.0f,%2.0f) microns, ecc = (%2.1f, %2.1f) degs', ...
        coneMosaicCenterPositionMM(1)*1e3, coneMosaicCenterPositionMM(2)*1e3, mosaicEccDegs(1), mosaicEccDegs(2)));
    
    ax = theAxesGrid{1,2};
    hold(ax, 'on');
    for k = 1:size(RGCRFPositionsMicrons,1)
        r = 0.5*RGCRFSpacingsMicrons(k);
        patch(ax,RGCRFPositionsMicrons(k,1)+r*xx, RGCRFPositionsMicrons(k,2)+r*yy,[0.8 0.8 0.8]);
    end
    axis(ax, 'equal'); axis(ax, 'square')
    set(ax, 'XLim', xRange, 'YLim', yRange, 'YTickLabel', {});
    xlabel(ax, 'microns');
    title(ax,sprintf('mRGC mosaic (cropped, with %2.0f micron padding for RF surrounds)', extraMicronsForSurroundCones));
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
