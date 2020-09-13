function theMidgetRGCmosaic = generateMRGCMosaicConnectedToConeMosaic(theConeMosaicMetaData, mRGCmosaicFile, mosaicParams, ...
    deconvolutionOpticsParams, visualizeSynthesizedParams, outputFile, exportsDir)
    
    % STEP 1. Retrieve regular hex cone mosaic metadata
    coneMosaicEccDegs = theConeMosaicMetaData.coneMosaicEccDegs;
    coneMosaicSizeMicrons = theConeMosaicMetaData.coneMosaicSizeMicrons;
    conePositionsMicrons = theConeMosaicMetaData.conePositionsMicrons;
    coneSpacingsMicrons = theConeMosaicMetaData.coneSpacingsMicrons;
    coneTypes = theConeMosaicMetaData.coneTypes;
    extraMicronsForSurroundCones = theConeMosaicMetaData.extraMicronsForSurroundCones;
    

    % STEP 2. Connect the cone mosaic patch to the centers of the midget RGC mosaic
    orphanRGCpolicy = mosaicParams.orphanRGCpolicy;
    maximizeConeSpecificity = mosaicParams.maximizeConeSpecificity;
    visualizeMosaicsToBeConnected = ~true;
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, midgetRGCconnectionMatrix] = ...
         connectMidgetRGCMosaicToConeMosaic(mRGCmosaicFile, mosaicParams.rgcMosaicPatchSizeMicrons, ...
         conePositionsMicrons, coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones, ...
         orphanRGCpolicy, maximizeConeSpecificity, visualizeMosaicsToBeConnected);

     
    % Visualize EXCLUSIVE connections to the RF centers. These are
    % determined solely by the relationship of the cone/mRGCRF densities
    visualizeRFcenterTiling = ~true;
    if (visualizeRFcenterTiling)
        subregionToVisualize.center = round(mosaicParams.rgcMosaicPatchEccMicrons);
        subregionToVisualize.size = coneMosaicSizeMicrons;
        visualizeCenterConnections(midgetRGCconnectionMatrix, RGCRFPositionsMicrons,...
                conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                coneMosaicEccDegs, subregionToVisualize, ...
                outputFile,exportsDir);
    end
    
    % STEP 3. Compute weighted connections to center/surround regions
    [midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
     synthesizedRFParams] = computeWeightedConeInputsToRGCCenterSurroundSubregions(...
            conePositionsMicrons,  coneTypes, ...
            midgetRGCconnectionMatrix, ...
            mosaicParams.rgcMosaicPatchEccMicrons, mosaicParams.rgcMosaicPatchSizeMicrons, ...
            deconvolutionOpticsParams, visualizeSynthesizedParams, exportsDir);
        
    % Visualize ALL (EXCLUSIVE+SHARED) connections to the RF centers.
    % The Shared cone connections are guided by weigthed connections, which
    % are determined by the CronerKaplan decolvolution model
    visualizeRFcenterTiling = ~true;
    if (visualizeRFcenterTiling)
        subregionToVisualize.center = round(mosaicParams.rgcMosaicPatchEccMicrons);
        subregionToVisualize.size = coneMosaicSizeMicrons;
        visualizeCenterConnections(midgetRGCconnectionMatrixCenter, RGCRFPositionsMicrons,...
                conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
                coneMosaicEccDegs, subregionToVisualize, ...
                outputFile,exportsDir);
    end
    
    
    % The midget RGC mosaic object (for now just the cone weights to the
    % center and surround regions)
    theMidgetRGCmosaic = struct(...
        'centerWeights', midgetRGCconnectionMatrixCenter, ...   % sparse matrix of weights for center cone signals, indexed according to the serialization order of the cone mosaic
        'surroundWeights', midgetRGCconnectionMatrixSurround, ...  % sparse matrix of weights for surround cone signals, indexed according to the serialization order of the cone mosaic
        'extraMicronsForSurroundCones', extraMicronsForSurroundCones, ...
        'synthesizedRFParams', synthesizedRFParams);
    
    visualizeRFs = ~true;
    if (visualizeRFs)
        % Visualize the generated retinal 2D RFs (video)
        plotlabOBJ = setupPlotLab();
        
        outputFile = sprintf('%s_RFexamples',outputFile);
        visualizeSubregions(1,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
            synthesizedRFParams.rfEccRadiusDegs, ...
            synthesizedRFParams.rfCenterPositionMicrons,  ...
            synthesizedRFParams.retinal.surroundCharacteristicRadiiMicrons,...
            conePositionsMicrons, coneSpacingsMicrons,  coneTypes, ...
            plotlabOBJ, outputFile, exportsDir);
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

    showConnectedConePolygon = ~true;
    visualizeRFs(coneMosaicEccDegs, zLevels, whichLevelsToContour, ...
             midgetRGCconnectionMatrix, RGCRFPositionsMicrons,...
             conePositionsMicrons, coneSpacingsMicrons, coneTypes, subregionToVisualize, ...
             displayEllipseInsteadOfContour, showConnectedConePolygon, plotlabOBJ,  outputFile, exportsDir);
         
end


function [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, midgetRGCconnectionMatrix] = ...
    connectMidgetRGCMosaicToConeMosaic(mRGCmosaicFile, rgcMosaicPatchSizeMicrons, conePositionsMicrons, ...
    coneSpacingsMicrons, coneTypes, extraMicronsForSurroundCones, orphanRGCpolicy, maximizeConeSpecificity, visualizedMosaics)

    % Load mRGC RF data
    load(mRGCmosaicFile, ...
        'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios');

    % Crop midget mosaic to the size and position of the cone mosaic
  	mRGCRFroi.center = 0.5*(min(conePositionsMicrons, [], 1) + max(conePositionsMicrons, [], 1));
    mRGCRFroi.size = max(conePositionsMicrons, [], 1) - min(conePositionsMicrons, [], 1);
    
    
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios] = ...
        cropRGCmosaic(RGCRFPositionsMicrons, RGCRFSpacingsMicrons,  desiredConesToRGCratios, mRGCRFroi);
    
    
    
    % Compute inputs to RGC RF centers
    visualizeConnectionProcess = ~true;
    [midgetRGCconnectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons] = computeConnectionMatrix(...
                RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, coneSpacingsMicrons, ...
                coneTypes, desiredConesToRGCratios, orphanRGCpolicy, maximizeConeSpecificity, ...
                visualizeConnectionProcess);
            
           
    % Compute RGC positions from connectivity here
    rgcsNum = size(midgetRGCconnectionMatrix,2);
    RGCRFPositionsMicronsFromConnectivity = zeros(rgcsNum,2);
    for iRGC = 1:size(RGCRFPositionsMicrons,1)
        weights = full(squeeze(midgetRGCconnectionMatrix(:, iRGC)));
        centerIndices = find(weights>0);
        RGCRFPositionsMicronsFromConnectivity(iRGC,:) = mean(conePositionsMicrons(centerIndices,:),1);
    end
    
    % Only keep RGCs within mRGCRFroi.center +/-
    % 0.5*rgcMosaicPatchSizeMicrons AND allowing space for surround cones
    minConePosition = min(conePositionsMicrons, [], 1);
    maxConePosition = max(conePositionsMicrons, [], 1);

    finalRGCindices = [];
    for rgcIndex = 1:size(RGCRFPositionsMicrons,1)
        distanceVector = abs(RGCRFPositionsMicronsFromConnectivity(rgcIndex,:) - mRGCRFroi.center);
        if ( ...
           (distanceVector(1) <= 0.5*rgcMosaicPatchSizeMicrons(1)) && (distanceVector(2) <= 0.5*rgcMosaicPatchSizeMicrons(2)) && ...
           (RGCRFPositionsMicronsFromConnectivity(rgcIndex,1)-extraMicronsForSurroundCones/2 >= minConePosition(1)) && ...
           (RGCRFPositionsMicronsFromConnectivity(rgcIndex,1)+extraMicronsForSurroundCones/2 <= maxConePosition(1)) && ...
           (RGCRFPositionsMicronsFromConnectivity(rgcIndex,2)-extraMicronsForSurroundCones/2 >= minConePosition(2)) && ...
           (RGCRFPositionsMicronsFromConnectivity(rgcIndex,2)+extraMicronsForSurroundCones/2 <= maxConePosition(2)) ...
           )
            finalRGCindices = cat(2, finalRGCindices, rgcIndex);
        end
    end
    
    midgetRGCconnectionMatrix = midgetRGCconnectionMatrix(:, finalRGCindices);
    RGCRFPositionsMicrons = RGCRFPositionsMicronsFromConnectivity(finalRGCindices,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(finalRGCindices);
    
    % Visualize connected mosaics
    if (visualizedMosaics)
        visualizeMosaicPatchesToBeConnected(conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            RGCRFPositionsMicrons, RGCRFSpacingsMicrons, mRGCRFroi.center)
    end
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