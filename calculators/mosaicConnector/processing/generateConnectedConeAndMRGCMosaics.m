function [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams)
    % STEP 1. Generate a regular hex cone mosaic patch with the desired eccentricity and size
    [theConeMosaic, coneMosaicEccDegs, coneMosaicSizeMicrons, conePositionsMicrons, coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones] = ...
        generateRegularHexMosaicPatch(...
            mosaicParams.rgcMosaicPatchEccMicrons, ...
            mosaicParams.rgcMosaicPatchSizeMicrons);
     
    % STEP 2. Connect the cone mosaic patch to the centers of the midget RGC mosaic
    orphanRGCpolicy = mosaicParams.orphanRGCpolicy;
    maximizeConeSpecificity = mosaicParams.maximizeConeSpecificity;
    visualizeMosaicsToBeConnected = ~true;
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, midgetRGCconnectionMatrix] = ...
         connectMidgetRGCMosaicToConeMosaic(mRGCmosaicFile, mosaicParams.rgcMosaicPatchSizeMicrons, ...
         conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
         orphanRGCpolicy, maximizeConeSpecificity, visualizeMosaicsToBeConnected);
     
    % Visualize connections to the RF centers
    visualizeRFcenterTiling = ~true;
    if (visualizeRFcenterTiling)
        subregionToVisualize.center = round(runParams.rgcMosaicPatchEccMicrons);
        subregionToVisualize.size = coneMosaicSizeMicrons;
        visualizeCenterConnections(midgetRGCconnectionMatrix, RGCRFPositionsMicrons,...
                conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                coneMosaicEccDegs, subregionToVisualize, ...
                runParams.outputFile,runParams.exportsDir);
    end
    
    % STEP 3. 
    [midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
     synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
            conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            RGCRFPositionsMicrons, midgetRGCconnectionMatrix, ...
            mosaicParams.rgcMosaicPatchEccMicrons, mosaicParams.rgcMosaicPatchSizeMicrons);
        
    % The midget RGC mosaic object (for now just the cone weights to the
    % center and surround regions)
    theMidgetRGCmosaic = struct(...
        'centerWeights', midgetRGCconnectionMatrixCenter, ...   % sparse matrix of weights for center cone signals, indexed according to the serialization order of the cone mosaic
        'surroundWeights', midgetRGCconnectionMatrixSurround, ...  % sparse matrix of weights for surround cone signals, indexed according to the serialization order of the cone mosaic
        'extraMicronsForSurroundCones', extraMicronsForSurroundCones);
    
    visualizeRFs = ~true;
    if (visualizeRFs)
        % Visualize the generated retinal 2D RFs (video)
        plotlabOBJ = setupPlotLab();

        outputFile = sprintf('%s_RFexamples',runParams.outputFile);
        visualizeSubregions(1,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
            synthesizedRFParams.rgcIndices,  synthesizedRFParams.eccDegs, ...
            synthesizedRFParams.centerPositionMicrons, synthesizedRFParams.retinal.centerRadiiDegs, ...
            synthesizedRFParams.retinal.surroundRadiiDegs,...
            conePositionsMicrons, coneSpacingsMicrons,  coneTypes, ...
            plotlabOBJ, outputFile, runParams.exportsDir);
    end
end


function visualizeCenterConnections(midgetRGCconnectionMatrix, RGCRFPositionsMicrons,...
            conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            coneMosaicEccDegs, subregionToVisualize, outputFile, exportsDir)
    %Visualize the connections to the RF centers
    zLevels = [0.3 1];
    whichLevelsToContour = [1];
    displayEllipseInsteadOfContour = false;
    
    figHeightInches = 15;
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
                'renderer', 'painters', ...
                'axesBox', 'on', ...
                'colorOrder', [0 0 0; 1 0 0.5], ...
                'axesTickLength', [0.015 0.01]/4,...
                'axesFontSize', 22, ...
                'figureWidthInches', figHeightInches/(subregionToVisualize.size(2))*(subregionToVisualize.size(1)), ...
                'figureHeightInches', figHeightInches);

    visualizeRFs(coneMosaicEccDegs, zLevels, whichLevelsToContour, ...
             midgetRGCconnectionMatrix, RGCRFPositionsMicrons,...
             conePositionsMicrons, coneSpacingsMicrons, coneTypes, subregionToVisualize, ...
             displayEllipseInsteadOfContour, plotlabOBJ,  outputFile, exportsDir);
         
end

function extraMicronsForSurroundCones = estimateMaxSurroundRadiusMicrons(eccentricityMicrons, sizeMicrons)
    w = WatsonRGCModel('generateAllFigures', false);
    posUnits = 'mm'; densityUnits = 'mm^2';
    coneEccMaxMicrons = max(abs([eccentricityMicrons+sizeMicrons/2; eccentricityMicrons-sizeMicrons/2]));
    
    % Cone spacing at the most eccentric position
    coneSpacingMicronsMax = 1e3 * w.coneRFSpacingAndDensityAtRetinalPositions(...
        coneEccMaxMicrons*1e-3, 'right', posUnits, densityUnits, ...
        'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
 
    % Compute RF params using the CronerKaplan model
    ck = CronerKaplanRGCModel('generateAllFigures', false, 'instantiatePlotLab', false);
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(coneSpacingMicronsMax, coneEccMaxMicrons);
    retinalSurroundRadiusDegsAt1overE = synthesizedRFParams.retinal.surroundRadiiDegs;
    extraDegsForSurroundCones = retinalSurroundRadiusDegsAt1overE*2;
    extraMicronsForSurroundCones = ceil(WatsonRGCModel.sizeDegsToSizeRetinalMicrons(extraDegsForSurroundCones, synthesizedRFParams.eccDegs));
end

function [theConeMosaic, coneMosaicEccDegs, coneMosaicSizeMicrons, conePositionsMicrons, coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones] = ...
    generateRegularHexMosaicPatch(eccentricityMicrons, sizeMicrons)

    % Estimate max RF surround radius
    extraMicronsForSurroundCones = estimateMaxSurroundRadiusMicrons(eccentricityMicrons, sizeMicrons);

    % Compute the cone mosaic FOV in degrees
    coneMosaicCenterPositionMM = eccentricityMicrons * 1e-3;
    coneMosaicSizeMicrons = sizeMicrons + 2*extraMicronsForSurroundCones*[1 1];
    coneMosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMM);
    fovDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(coneMosaicSizeMicrons, sqrt(sum((coneMosaicCenterPositionMM*1e3).^2,2.0)));
    
    % Determine the median cone spacing with the patch
    whichEye = 'right';
    coneSpacingMicrons = medianConeSpacingInPatch(whichEye, eccentricityMicrons, coneMosaicSizeMicrons);
    
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
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
    % Cone spacings: all the same
    coneSpacingsMicrons = ones(size(conePositionsMicrons,1),1) * coneSpacingMicrons;
    % Cone types
    coneTypes = cmStruct.coneTypes;
 
end

function [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, midgetRGCconnectionMatrix] = ...
    connectMidgetRGCMosaicToConeMosaic(mRGCmosaicFile, rgcMosaicPatchSizeMicrons, conePositionsMicrons, ...
    coneSpacingsMicrons, coneTypes, orphanRGCpolicy, maximizeConeSpecificity, visualizedMosaics)

    % Load mRGC RF data
    load(mRGCmosaicFile, ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');

    % Crop midget mosaic to the size and position of the cone mosaic, leaving enough space for the surround cones
  	mRGCRFroi.center = 0.5*(min(conePositionsMicrons, [], 1) + max(conePositionsMicrons, [], 1));
    mRGCRFroi.size = max(conePositionsMicrons, [], 1) - min(conePositionsMicrons, [], 1);
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
        cropRGCmosaic(RGCRFPositionsMicrons, RGCRFSpacingsMicrons,  desiredConesToRGCratios, mRGCRFroi);
    
    % Visualize mosaics to be connected
    if (visualizedMosaics)
        visualizeMosaicPatchesToBeConnected(conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            RGCRFPositionsMicrons, RGCRFSpacingsMicrons, mRGCRFroi.center)
    end
    
    % Compute inputs to RGC RF centers
    visualizeConnectionProcess = ~true;
    [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = computeConnectionMatrix(...
                RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
                coneTypes, desiredConesToRGCratios, orphanRGCpolicy, maximizeConeSpecificity, ...
                visualizeConnectionProcess);
            
    % Only keep RGCs within mRGCRFroi.center +/- 0.5*rgcMosaicPatchSizeMicrons
    finalRGCindices = [];
    for rgcIndex = 1:size(RGCRFPositionsMicrons,1)
        distanceVector = abs(RGCRFPositionsMicrons(rgcIndex,:) - mRGCRFroi.center);
        if (distanceVector(1) <= 0.5*rgcMosaicPatchSizeMicrons(1)) && (distanceVector(2) <= 0.5*rgcMosaicPatchSizeMicrons(2))
            finalRGCindices = cat(2, finalRGCindices, rgcIndex);
        end
    end
    midgetRGCconnectionMatrix = midgetRGCconnectionMatrix(:, finalRGCindices);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(finalRGCindices,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(finalRGCindices);
end

function visualizeMosaicPatchesToBeConnected(conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons, coneMosaicCenterPositionMicrons)
    
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
    mosaicEccDegs = WatsonRGCModel.rhoMMsToDegs(coneMosaicCenterPositionMicrons*1e-3);
    title(ax,sprintf('reg hex cone mosaic patch (x,y) = (%2.0f,%2.0f) microns, ecc = (%2.1f, %2.1f) degs', ...
        coneMosaicCenterPositionMicrons(1), coneMosaicCenterPositionMicrons(2), mosaicEccDegs(1), mosaicEccDegs(2)));
    
    ax = theAxesGrid{1,2};
    hold(ax, 'on');
    for k = 1:size(RGCRFPositionsMicrons,1)
        r = 0.5*RGCRFSpacingsMicrons(k);
        patch(ax,RGCRFPositionsMicrons(k,1)+r*xx, RGCRFPositionsMicrons(k,2)+r*yy,[0.8 0.8 0.8]);
    end
    axis(ax, 'equal'); axis(ax, 'square')
    set(ax, 'XLim', xRange, 'YLim', yRange, 'YTickLabel', {});
    xlabel(ax, 'microns');
    title(ax,sprintf('mRGC mosaic'));
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