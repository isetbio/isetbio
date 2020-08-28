function runPhaseX(runParams)
    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Load/Recompute connected mosaics and the optics
    recomputeConeMosaic = ~true;
    recomputeOptics = ~true;
    
    % mRGC mosaic: whether to re-generate it
    recomputeRGCmosaic = true;
    % mRGC mosaic: whether to visualize the synthesized RF params
    visualizeSynthesizedParams = true;
    
    % Compute cone mosaic responses
    recomputeConeMosaicResponses = ~true;
    recomputeNullResponses = ~true;
    
    % Responses directory
    saveDir = runParams.responseFilesDir;
    if (~isfolder(saveDir)) 
        if ((~recomputeConeMosaicResponses)&&(~recomputeNullResponses))
            % This directory should exist if we are computing RGC responses
            error('Could not find cone mosaic responses directory ''%s''.\n', saveDir);
        end
        % We are computing cone responses, so we may need to generate the
        % responses directory. Ask permission to create it.
        fprintf('Will generate responses saveDir : ''%s''\n', saveDir);
        fprintf('Hit enter to proceed: ');
        pause;
        mkdir(saveDir);
    end
    
    
    % Figure exports dir
    figExportsDir = runParams.exportsDir;
    if (~isfolder(figExportsDir))
        % Fig exports directory does not exist. Ask permission to
        % create it.
        fprintf('Will generate figExportsDir: ''%s''\n', figExportsDir);
        fprintf('Hit enter to proceed: ');
        pause;
        mkdir(figExportsDir);
    end

    
    % Generate the mosaics and the optics
    [theConeMosaic, theMidgetRGCmosaic, theOptics, opticsPostFix, PolansSubjectID] = ...
        mosaicsAndOpticsForEccentricity(runParams, recomputeConeMosaic, ...
        recomputeRGCmosaic, recomputeOptics, saveDir, ...
        visualizeSynthesizedParams);

    % Display PSFs
    displayPSFs = ~true;
    if (displayPSFs)
        eccDegs = WatsonRGCModel.rhoMMsToDegs(1e-3*runParams.rgcMosaicPatchEccMicrons);
        visualizePSFs(theOptics, eccDegs(1), eccDegs(2));
    end
    
    
    % Visual stimulation parameters
    LMScontrast = [0.1 0.1 0.0];
    minSF = 0.1;
    maxSF = 100;
    sfsNum = 15;
    spatialFrequenciesCPD = logspace(log10(minSF), log10(maxSF),sfsNum);
    
    stimulusFOVdegs = 2.0;
    minPixelsPerCycle = 8;
    stimulusPixelsNum = stimulusFOVdegs*minPixelsPerCycle;
    temporalFrequency = 4.0;
    stimDurationSeconds = 0.5;
    instancesNum = 64;
    
    % Visualized cells
    targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves = [1:4];
    
    stimColor = struct(...
        'backgroundChroma', [0.3, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', LMScontrast);
    
    stimTemporalParams = struct(...
        'temporalFrequencyHz', temporalFrequency, ...
        'stimDurationSeconds', stimDurationSeconds);
    
    stimSpatialParams = struct(...
        'type', 'driftingGrating', ...
        'fovDegs', stimulusFOVdegs,...
        'pixelsNumPerSF', stimulusPixelsNum, ...
        'minPixelsNum', 1024, ...
        'gaborPosDegs', [0 0], ...
        'gaborSpatialFrequencyCPD', 0, ...
        'gaborSigmaDegs', Inf, ... %stimulusFOVdegs/(2*4), ...%Inf, ...
        'gaborOrientationDegs', 0, ...
        'deltaPhaseDegs', []);
    
    % Signal to the RGCs
    rgcInputSignal = 'isomerizations';
    %rgcInputSignal = 'photocurrents';
    
    if (recomputeConeMosaicResponses)
        computeConeResponses(runParams, ...
            stimColor,  stimTemporalParams, stimSpatialParams, ...
            theConeMosaic, theOptics, ...
            recomputeNullResponses, ...
            instancesNum, ...
            spatialFrequenciesCPD, ...
            opticsPostFix, PolansSubjectID, ...
            saveDir, ...
            'saveCornealStimulusSequence', ~true, ...
            'saveRetinalStimulusSequence', ~true);
    else
        visualizeAllSpatialFrequencyTuningCurves = true;
        visualizeResponseComponents = ~true;
        visualizeRetinalContrasts = ~true;
        coVisualizeRetinalStimulusWithMosaics = ~true;
        visualizeMeanConeMosaicResponseAsAMovie = false;
        visualizeRGCTemporalResponsesAtRGCPositions = ~true;
        visualizeRGCSFTuningsAtRGCPositions = true;
        visualizePatchStatistics = true;
        
        computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
            rgcInputSignal, spatialFrequenciesCPD, LMScontrast, ...
            stimSpatialParams, stimTemporalParams, ...
            saveDir, figExportsDir, ...
            visualizeRGCTemporalResponsesAtRGCPositions, visualizeRGCSFTuningsAtRGCPositions, ...
            visualizeAllSpatialFrequencyTuningCurves, visualizeResponseComponents, ...
            visualizeRetinalContrasts, visualizeMeanConeMosaicResponseAsAMovie, ...
            targetRGCsForWhichToVisualizeSpatialFrequencyTuningCurves, visualizePatchStatistics, ...
            'coVisualizeRetinalStimulusWithMosaics', coVisualizeRetinalStimulusWithMosaics, ...
            'coVisualizedRetinalStimulusSpatialFrequency', 35, ...
            'coVisualizedRetinalStimulusConeContrast', LCONE_ID);
            
    end
end