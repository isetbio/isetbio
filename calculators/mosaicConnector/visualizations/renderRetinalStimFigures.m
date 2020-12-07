function coVisualizedRetinalStimulus = renderRetinalStimFigures(mFile, spatialFrequenciesCPD, ...
    coVisualizeRetinalStimulusWithMosaics, coVisualizedRetinalStimulusData, ...
    visualizeRetinalContrasts, visualizeMeanConeMosaicResponseAsAMovie, ...
    saveDir, runParams, LMScontrast, opticsPostFix, PolansSubjectID, ...
    figExportsDir)

    global LCONE_ID
    global MCONE_ID
    global SCONE_ID
    
    % Flags indicating whether the file containst corneal/retinal stimulus
    % sequence info
    cornealStimulusDataIsAvailable = mFile.saveCornealStimulusSequence;
    retinalStimulusDataIsAvailable = mFile.saveRetinalStimulusSequence;
    
    % Retrieve the time axes
    responseTimeAxis = mFile.responseTimeAxis;
    stimulusTimeAxis = mFile.stimulusTimeAxis;
            
    % Retrieve background retinal LMS excitations
    if (retinalStimulusDataIsAvailable)
        theRetinalLMSexcitationsSequenceNull = mFile.theRetinalLMSexcitationsSequence;
        backgroundLMSexcitations = mean(mean(mean(theRetinalLMSexcitationsSequenceNull,1),2),3);
    end
    
    % Stimulus to be co-visualized under the cone/RGC mosaic in target RGCs
    coVisualizedRetinalStimulus = [];
    
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        gaborSpatialFrequencyCPD = spatialFrequenciesCPD(sfIndex);
        
        % Assemble the test response filename
        theTestResponseFileName = sprintf('%s_%2.1fCPD.mat',testResponseFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID), gaborSpatialFrequencyCPD);
        
        % Open data file for read-only
        mFile = matfile(fullfile(saveDir, theTestResponseFileName), 'Writable', false);
        
        if (coVisualizeRetinalStimulusWithMosaics) && (retinalStimulusDataIsAvailable)
            switch (coVisualizedRetinalStimulusData.coneContrast)
                case LCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 1;  % visualize the L-cone contrast
                case MCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 2;  % visualize the M-cone contrast
                case SCONE_ID
                        visualizedRetinalStimulusConeContrastIndex = 3;  % visualize the S-cone contrast
            end
            [~,visualizedRetinalSFindex] = min(abs(spatialFrequenciesCPD-coVisualizedRetinalStimulusData.spatialFrequency));
          
            if (sfIndex == visualizedRetinalSFindex)    
                theRetinalLMSexcitations = mFile.theRetinalLMSexcitationsSequence(coVisualizedRetinalStimulusData.frameIndex,:,:,:);
                theRetinalLMScontrast = bsxfun(@times, bsxfun(@minus, theRetinalLMSexcitations, backgroundLMSexcitations), 1./backgroundLMSexcitations);
                coVisualizedRetinalStimulus.retinalContrastImage = squeeze(theRetinalLMScontrast(1,:,:,visualizedRetinalStimulusConeContrastIndex));
                coVisualizedRetinalStimulus.spatialSupportMicrons = mFile.retinalSpatialSupportMicrons;
            end
        end
        
        
        if (visualizeRetinalContrasts && retinalStimulusDataIsAvailable)
            % The test retinal LMS excitations 
            theRetinalLMSexcitationsSequence = mFile.theRetinalLMSexcitationsSequence;

            % Compute retinal contrast
            theRetinalLMScontrastSequence = bsxfun(@minus, theRetinalLMSexcitationsSequence, backgroundLMSexcitations);
            theRetinalLMScontrastSequence = bsxfun(@times, theRetinalLMScontrastSequence, 1./backgroundLMSexcitations);

            % Extract retinal L and M contrast over central region
            colsNum = size(theRetinalLMScontrastSequence,3);
            margin = round(colsNum/5);
            colsToUse = margin:2:(colsNum-margin);
            colsToUseNum = numel(colsToUse);
            midRow = round(0.5*size(theRetinalLMScontrastSequence,2));
            
            if (sfIndex == 1)
                roiLcontrast = cell(1, numel(spatialFrequenciesCPD));
                roiMcontrast = cell(1, numel(spatialFrequenciesCPD));
            end
            
            roiLcontrast{sfIndex} = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,1)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);
            roiMcontrast{sfIndex} = reshape(squeeze(theRetinalLMScontrastSequence(:,midRow,colsToUse,2)), [1 size(theRetinalLMScontrastSequence,1)*colsToUseNum]);

            visualizeRetinalContrastProfiles(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast, figExportsDir);
        end
        
        % Display the mean presynaptic response
        if (visualizeMeanConeMosaicResponseAsAMovie && cornealStimulusDataIsAvailable)
            % Retrieve the stimulus RGB sequence
            theStimulusRGBsequence = mFile.theStimulusSRGBsequence;
            theMeanPresynapticResponses = squeeze(mean(thePresynapticResponses,1));
            visualizeStimResponseMovie(responseTimeAxis, stimulusTimeAxis, stimSpatialParams, theStimulusRGBsequence, ...
                theMeanPresynapticResponses, theConeMosaic);
        end
    end % sfIndex
    
    if (visualizeRetinalContrasts) && (retinalStimulusDataIsAvailable)
        visualizeRetinalLMcontrastCorrelation(spatialFrequenciesCPD, roiLcontrast, roiMcontrast, ...
           LMScontrast, opticsPostFix, PolansSubjectID, true, figExportsDir);
    end
    
end