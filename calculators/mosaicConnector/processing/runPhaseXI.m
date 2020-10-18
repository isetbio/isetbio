function runPhaseXI(runParams)

    % Figure exports dir
    figExportsDir = runParams.exportsDir;
    
    stimulusType = 'drifting gratings';
    switch stimulusType
        case 'drifting gratings'
            plotChromaticTuning(runParams, figExportsDir );
    end
end


function plotChromaticTuning(runParams, figExportsDir)

    % Responses directory
    saveDir = runParams.responseFilesDir;
    
    [~, ~, opticsPostFix, PolansSubjectID] = mosaicsAndOpticsFileName(runParams);
    
    % Import the computed spatial transfer functions
    LMangles = [-45 0 45 90];
    
    for iAngle = 1:numel(LMangles)
        switch LMangles(iAngle)
            case -45
                LMScontrast = [0.1 -0.1 0.0];
            case 0 
                LMScontrast = [0.1 0.0 0.0];
            case 45
                LMScontrast = [0.1 0.1 0.0];
            case 90
                LMScontrast = [0.0 0.1 0.0];
        end
        
        theSpatialTransferFunctionsFileName = fullfile(saveDir, spatialTransferFunctionsFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID));
        if (LMScontrast(1) == LMScontrast(2))
            load(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude', ...
             'theConeMosaic', 'theMidgetRGCmosaic');
        else
            load(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude');
        end
        
        % Scale tuning by the RMS cone contrast
        LMtuning(iAngle, :,:) = responseAmplitude / (norm(LMScontrast)*100);
    end % iAngle
              
    visualizedSet = 'data and fit';
    %visualizedSet = 'fit and isoresponse';
    
    [~,rgcsNum, sfsNum] = size(LMtuning);
    
    % Which SF to use for the chromatic tuning
    targetSFindex = 1;
    
    for iRGC = 1:rgcsNum
        
        [xFullData, yFullData, xCosineFit, yCosineFit, xIsoResponseFit, yIsoResponseFit, phaseDegs, peakResponse, baserate] = ...
            fitCosineToLMplaneResponses(squeeze(LMtuning(:,iRGC, targetSFindex)), LMangles);
        
        switch visualizedSet
            case 'data and fit'
                LMplane.x(iRGC,:) = xFullData;
                LMplane.y(iRGC,:) = yFullData;
                LMplane.xFit(iRGC,:) = xCosineFit;
                LMplane.yFit(iRGC,:) = yCosineFit;
            
            case 'fit and isoresponse'
                LMplane.x(iRGC,:) = xCosineFit;
                LMplane.y(iRGC,:) = yCosineFit;
                LMplane.xFit(iRGC,:) = xIsoResponseFit;
                LMplane.yFit(iRGC,:) = yIsoResponseFit;
                
            otherwise
                error('Unknown data set');
        end
    end
    
    
    plotXaxisScaling = 'linear';
    plotType = 'LMplaneTuning';
    

    maxSpikeRateModulation = 1.0*max([prctile(LMplane.x(:),99) prctile(LMplane.y(:),99)])

    spikeRateTicks = 2:2:20;
    coVisualizedRetinalStimulus = [];
    exportFig = true;
    figureName = 'LMplaneTuningAll';
    
    visualizeRGCmosaicWithResponses(1000, theConeMosaic, plotXaxisScaling, plotType, ...
                        LMplane.x,LMplane.y, ...
                        LMplane.xFit,  LMplane.yFit, ...
                        [], [], ...
                        runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                        theMidgetRGCmosaic,  'centers', maxSpikeRateModulation, spikeRateTicks, ...
                        coVisualizedRetinalStimulus, ...
                        figureName, [nan nan nan], opticsPostFix, PolansSubjectID, ...
                        [], false, ...
                        exportFig, figExportsDir);
                    
    for targetRGC = 1:rgcsNum
        figureName = sprintf('LMplaneTuning%d', targetRGC);

        visualizeRGCmosaicWithResponses(1000+targetRGC, theConeMosaic, plotXaxisScaling, plotType, ...
                        LMplane.x,LMplane.y, ...
                        LMplane.xFit,  LMplane.yFit, ...
                        [], [], ...
                        runParams.rgcMosaicPatchEccMicrons, runParams.rgcMosaicPatchSizeMicrons, ...
                        theMidgetRGCmosaic,  'centers', maxSpikeRateModulation, spikeRateTicks, ...
                        coVisualizedRetinalStimulus, ...
                        figureName, [nan nan nan], opticsPostFix, PolansSubjectID, ...
                        targetRGC, false, ...
                        exportFig, figExportsDir);
    end
end
