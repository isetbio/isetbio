function runPhaseX(runParams)

    photocurrentResponseDataFile = 'photocurrentResponses2.mat';
    
    recomputePhotocurrents = ~true;
    if (recomputePhotocurrents)
        instancesNum = 2;
        stimDurationSeconds = 0.1;
        computePhotocurrents(runParams, instancesNum, stimDurationSeconds, photocurrentResponseDataFile);
    else
        computeRGCresponses(runParams, photocurrentResponseDataFile);
    end
end

function computeRGCresponses(runParams, photocurrentResponseDataFile)

    mFile = matfile(photocurrentResponseDataFile, 'Writable', false);
    fprintf('\nImporting the cone mosaic ...');
    theConeMosaic = mFile.theConeMosaic;
    fprintf('Done !\n');
    fprintf('\nImporting the RGC mosaic ...');
    theMidgetRGCmosaic = mFile.theMidgetRGCmosaic;
    fprintf('Done !\n');
    fprintf('\nVisualizing the RGC mosaic ...');
    visualizeConeAndRGCmosaics(theConeMosaic, runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, theMidgetRGCmosaic);
    fprintf('Done !\n');
    pause
    
    theOISequence = mFile.theOIsequence;
    thePhotocurrents = mFile.photocurrents;
    
    size(photocurrents)
    class(photocurrents)
end

function visualizeConeAndRGCmosaics(theConeMosaic, eccentricityMicrons, sizeMicrons, theMidgetRGCmosaic)
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the mosaic center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, eccentricityMicrons);
   
    % Cone types
    coneTypes = cmStruct.coneTypes;
    coneDiameterMicrons = cmStruct.coneApertures * theConeMosaic.micronsPerDegree;
    coneSpacingsMicrons = 1.0/0.7 * coneDiameterMicrons;
    
    conesNum = size(theMidgetRGCmosaic.centerWeights,1);
    rgcsNum = size(theMidgetRGCmosaic.centerWeights,2);
    
    rgcPos = zeros(rgcsNum,2);

    % Sampling for contours
    deltaX = 0.25;
    xAxis = (eccentricityMicrons(1)-sizeMicrons(1)/2): deltaX: (eccentricityMicrons(1)+sizeMicrons(1)/2);
    yAxis = (eccentricityMicrons(2)-sizeMicrons(2)/2): deltaX: (eccentricityMicrons(2)+sizeMicrons(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    hFig = figure(99); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.05, ...
            'bottomMargin', 0.05, ...
            'rightMargin', 0.03, ...
            'topMargin', 0.1);
        
    theAxesGrid = theAxesGrid{1,1};
    xLims = [xAxis(1) xAxis(end)];
    yLims = [yAxis(1) yAxis(end)];
    ylabel(theAxesGrid, 'microns');
    set(theAxesGrid, 'CLim', [0 1], 'XLim', xLims, 'YLim', yLims);
    hold(theAxesGrid, 'on');
    colormap(theAxesGrid, brewermap(512, 'greys'));
    
    % Display cones
    LconeIndices = find(coneTypes == LCONE_ID);
    MconeIndices = find(coneTypes == MCONE_ID);
    SconeIndices = find(coneTypes == SCONE_ID);
    scatter(theAxesGrid,conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    scatter(theAxesGrid,conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0.5 0.9 0.5]);
    scatter(theAxesGrid,conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1.0]);
    
    
    for mRGCindex = 1:rgcsNum
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, mRGCindex)));
        centerIndices = find(centerWeights>0);
        
        % Use binary weights for visualization
        centerWeightsForVisualization = centerWeights;
        centerWeightsForVisualization(centerIndices) = 1;
        
        surroundWeights = full(theMidgetRGCmosaic.surroundWeights(:, mRGCindex));
        surroundIndices = find(surroundWeights>0);
        
        % Generate RF centers of RGCs based on cone positions and connection matrix
        theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(...
            centerWeightsForVisualization, conePositionsMicrons, coneSpacingsMicrons, X,Y);
        
        if (isempty(theRF))
            fprintf(2,'No cone inputs to this RF -> No visualization\n');
            continue;
        end
        
        zLevels = [0.3 1];
        whichLevelsToContour = 1;
        fitEllipse = false;
        
        C = contourc(xAxis, yAxis, theRF, zLevels);
        fillRFoutline(theAxesGrid, C, zLevels, whichLevelsToContour, fitEllipse);
        displayConnectedConesPolygon(theAxesGrid, centerIndices, conePositionsMicrons);
        
        drawnow;
        
        fprintf('mRGC %d has %d inputs in its center and %d inputs in its surround\n', ...
            mRGCindex, numel(centerIndices), numel(surroundIndices));
        
        rgcPos(mRGCindex,:) = mean(conePositionsMicrons(centerIndices,:),1);
        
    end % mRGCindex
    
    figure(100); clf;
    plot(rgcPos(:,1), rgcPos(:,2), 'ko'); hold on
end


function computePhotocurrents(runParams, instancesNum, stimDurationSeconds, photocurrentResponseDataFile)
    % Params for the connected mRGC mosaic
    mosaicParams.rgcMosaicPatchEccMicrons = runParams.rgcMosaicPatchEccMicrons;
    mosaicParams.rgcMosaicPatchSizeMicrons = runParams.rgcMosaicPatchSizeMicrons;
    mosaicParams.orphanRGCpolicy = runParams.orphanRGCpolicy;
    mosaicParams.maximizeConeSpecificity = runParams.maximizeConeSpecificity;
    
    % Location of file with the full (50 x 50 deg) ecc-based mRGCRF lattice
    mRGCmosaicFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
    
    % Generate cone mosaic and connected mRGC mosaic patches
    [theConeMosaic, theMidgetRGCmosaic] = generateConnectedConeAndMRGCMosaics(mRGCmosaicFile, mosaicParams);
    
    % Generate optics appropriate for the RGC mosaic eccentricity
    wavelengthSampling = theConeMosaic.pigment.wave;
    pupilDiameterMM = 3.0;
    wavelengthsListToCompute = wavelengthSampling;
    wavefrontSpatialSamples = 1001;
    micronsPerDegree = []; % empty so as to compute for each eccentricity
    imposedRefractionErrorDiopters = 0;
    deltaEcc = 1;
    eccXrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(1))*[1 1];
    eccYrangeDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*mosaicParams.rgcMosaicPatchEccMicrons(2))*[1 1];
    PolansSubjectID = 9;
    
    [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, wavefrontSpatialSamples, ...
                eccXrangeDegs, eccYrangeDegs, deltaEcc);

    theOI = theOIs{1,1,1};
    visualizePSFs(theOI, eccXrangeDegs(1), eccYrangeDegs(1));
    
    % Generate sine-wave scene as realized on a display
    stimColor = struct(...
        'backgroundChroma', [0.31, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', [0.1 0.0 0.0]);
    
    stimTemporalParams = struct(...
        'temporalFrequencyHz', 10.0, ...
        'stimDurationSeconds', stimDurationSeconds);
    
    stimSpatialParams = struct(...
        'fovDegs', 2.0,...
        'pixelsNum', 256, ...
        'gaborPosDegs', [0 0], ...
        'gaborSpatialFrequencyCPD', 4.0, ...
        'gaborSigma', Inf, ...
        'gaborOrientationDegs', 0, ...
        'deltaPhaseDegs', 30);
    
    % Generate scenes corresponding to each spatial phase of a drifting grating
    theSceneFrames = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling);
    
    % Generate the OISequence
    fprintf('\nComputing the oiSequence ...');
    tic
    theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOI, stimSpatialParams, stimTemporalParams);
    %theOIsequence.visualize('montage', 'backendrenderer', 'figure');
    fprintf('Done in %2.1f minutes!\n', toc/60);
    
    % Generate eye movements for the oiSequence
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    emPaths = zeros(instancesNum, eyeMovementsNum, 2);
        
    % Compute isomerizations and photocurrents
    fprintf('\nComputing the mosaic response ...');
    tic
    [isomerizations, photocurrents, osLinearFilters] = ...
        theConeMosaic.computeForOISequence(theOIsequence, ...
        'emPaths', emPaths, 'currentFlag', true, ...
        'workerID', 1);
    fprintf('Done in %2.1f minutes!\n', toc/60);
    
    fprintf('\nExporting data ...');
    tic
    save(photocurrentResponseDataFile, ...
        'theConeMosaic', ...
        'theMidgetRGCmosaic', ...
        'theOIsequence', ...
        'isomerizations', ...
        'photocurrents', ...
        '-v7.3');
    fprintf('Done in %2.1f minutes!\n', toc/60);
end

