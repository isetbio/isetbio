function runPhaseXI(runParams)

    % Figure exports dir
    figExportsDir = runParams.exportsDir;
    
    stimulusType = 'drifting gratings';
    switch stimulusType
        case 'drifting gratings'
            visualizedSet = 'data and fit';
            plotChromaticTuning(runParams, figExportsDir,visualizedSet );
            visualizedSet = 'fit and isoresponse';
            plotChromaticTuning(runParams, figExportsDir, visualizedSet );
    end
end


function plotChromaticTuning(runParams, figExportsDir, visualizedSet)

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
             
    
    [~,rgcsNum, sfsNum] = size(LMtuning);
    
    % Which SF to use for the chromatic tuning
    targetSFindex = 1;
    
    % isoresponse thresholds to attain 50% of max response
    isoResponseLevelPercentOfMax = 0.2;
    
    
    for iRGC = 1:rgcsNum
        
        [xFullData, yFullData, xCosineFit, yCosineFit, xIsoResponseFit, yIsoResponseFit, phaseDegs(iRGC), peakResponse, baserate] = ...
            fitCosineToLMplaneResponses(squeeze(LMtuning(:,iRGC, targetSFindex)), LMangles, isoResponseLevelPercentOfMax);
        
        switch visualizedSet
            case 'data and fit'
                LMplane.x(iRGC,:) = xFullData;
                LMplane.y(iRGC,:) = yFullData;
                LMplane.xFit(iRGC,:) = xCosineFit;
                LMplane.yFit(iRGC,:) = yCosineFit;
                figureName = 'LMplaneTuningAll_Response';
                
            case 'fit and isoresponse'
                LMplane.x(iRGC,:) = xCosineFit;
                LMplane.y(iRGC,:) = yCosineFit;
                LMplane.xFit(iRGC,:) = xIsoResponseFit;
                LMplane.yFit(iRGC,:) = yIsoResponseFit;
                figureName = 'LMplaneTuningAll_IsoResponseContrast';
                
            otherwise
                error('Unknown data set');
        end
    end
    
    hFig = figure(100);
    set(hFig, 'Color', [1 1 1]);
    phaseDegs = mod(phaseDegs+360,180);
    histogram(phaseDegs,0:5:180);
    axis 'square'
    set(gca, 'XTick', 0:30:180, 'XLim', [0 180]);
    grid on;
    xlabel('LM angle');
    ylabel('count');
    set(gca, 'FontSize', 20);
    pause
    
    plotXaxisScaling = 'linear';
    plotType = 'LMplaneTuning';
    

    maxSpikeRateModulation = 1.0*max([prctile(LMplane.x(:),99) prctile(LMplane.y(:),99)])

    spikeRateTicks = 2:2:20;
    coVisualizedRetinalStimulus = [];
    exportFig = true;
    
    % All RGCs
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
                    
    % Target RGCs
    
    for iRGC = 1:numel(runParams.targetRGCsForWhichToVisualizeTuningCurves)
        
        targetRGC = runParams.targetRGCsForWhichToVisualizeTuningCurves(iRGC);
        
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
