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
    
    % Signal to the RGCs
    rgcInputSignal = 'isomerizations';
    %rgcInputSignal = 'photocurrents';
    
    % Visualized cells
    targetRGCsForWhichToVisualizeTuningCurves =  [1:4]; %[20 21 22 23 24]; %[71 80 63 116]; 
   
    % Chromatic stimulus params
    LMScontrast = [0.1 0.1 0.0];
    stimColor = struct(...
        'backgroundChroma', [0.3, 0.31], ...
        'meanLuminanceCdPerM2', 40, ...
        'lmsContrast', LMScontrast);
    
    % Generate temporal and spatial stimulus params
    stimulusType = 'drifting gratings';
    instancesNum = 64;
    [stimTemporalParams, stimSpatialParams] = generateSpatioTemporalStimulationParams(stimulusType);
    
    % Visualization options
    visualizeAllTuningCurves = true;
    visualizeResponseComponents = true;
    visualizeRetinalContrasts = ~true;
    coVisualizeRetinalStimulusWithMosaics = ~true;
    visualizeMeanConeMosaicResponseAsAMovie = false;
    visualizeRGCTemporalResponsesAtRGCPositions = ~true;
    visualizeTuningAtRGCPositions = true;
    visualizePatchStatistics = true;
        
    if (recomputeConeMosaicResponses)
        switch stimulusType
            case 'drifting gratings'
            computeConeResponsesToDriftingGratings(runParams, ...
                stimColor,  stimTemporalParams, stimSpatialParams, ...
                theConeMosaic, theOptics, ...
                recomputeNullResponses, ...
                instancesNum, ...
                opticsPostFix, PolansSubjectID, saveDir, ...
                'saveCornealStimulusSequence', ~true, ...
                'saveRetinalStimulusSequence', ~true);
        end
        
    else  
        
        % Using Weber response representation ensures that the cone
        % responses are proportional to cone Weber contrast. 
        % That way a 10% L+M cone contrast will active L-cone center RGCs
        % to the same extent as it will activate M-cone center RGCs. If we
        % do not use this the relative response of L-cone center to M-cone center RGCs
        % depends on the mean cone response to the background.
        % Phototransductios does this transformation so that
        % pulses of same contrast have the same photocurrent modulation
        % independent on the background activation. But since we are
        % working on the cone isomerizations signal here, we have to do
        % this contrast transformation to have RGC activations that are
        % proportional to stimuls contrast.
        useWeberResponseRepresentation = true;
        
        switch stimulusType
            case 'drifting gratings'
            computeRGCresponsesToDriftingGratings(runParams, theConeMosaic, theMidgetRGCmosaic, ...
                rgcInputSignal, LMScontrast, ...
                stimSpatialParams, stimTemporalParams, ...
                useWeberResponseRepresentation, ...
                saveDir, figExportsDir, ...
                visualizeRGCTemporalResponsesAtRGCPositions, visualizeTuningAtRGCPositions, ...
                visualizeAllTuningCurves, visualizeResponseComponents, ...
                visualizeRetinalContrasts, visualizeMeanConeMosaicResponseAsAMovie, ...
                targetRGCsForWhichToVisualizeTuningCurves, visualizePatchStatistics, ...
                'coVisualizeRetinalStimulusWithMosaics', coVisualizeRetinalStimulusWithMosaics, ...
                'coVisualizedRetinalStimulusSpatialFrequencyCPD', 35, ...
                'coVisualizedRetinalStimulusConeContrast', LCONE_ID);
        end   
    end
end

