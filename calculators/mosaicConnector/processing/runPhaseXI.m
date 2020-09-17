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
    LMangles = [0 45 90 135];
    
    for iAngle = 1:numel(LMangles)
        switch LMangles(iAngle)
            case 0 
                LMScontrast = [0.1 0.0 0.0];
            case 45
                LMScontrast = [0.1 0.1 0.0];
            case 90
                LMScontrast = [0.0 0.1 0.0];
            case 135
                LMScontrast = [-0.1 0.1 0.0];
        end
        rmsContrast = sqrt(sum(LMScontrast.^2));
        theSpatialTransferFunctionsFileName = fullfile(saveDir, spatialTransferFunctionsFilename(runParams, LMScontrast, opticsPostFix, PolansSubjectID));
        if (LMScontrast(1) == LMScontrast(2))
            load(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude', ...
             'theConeMosaic', 'theMidgetRGCmosaic');
        else
            load(theSpatialTransferFunctionsFileName, 'spatialFrequenciesCPD', 'responseAmplitude');
        end
        
        % Scale tuning by the RMS contrast
        LMtuning(iAngle, :,:) = responseAmplitude / (rmsContrast*100);
    end % iAngle
              
    targetSFindex = 1;
    for iRGC = 1:size(LMtuning,2)
        x = squeeze(LMtuning(:,iRGC, targetSFindex)) .* cosd(LMangles(:));
        y = squeeze(LMtuning(:,iRGC, targetSFindex)) .* sind(LMangles(:));
        % make it symmetric
        x = cat(1, x, -x);
        y = cat(1, y, -y);
        [xEllipse,  yEllipse, semiAxes, rfCenter, noFit] = fitEllipseToContour(x,y);
        if (noFit)
            xEllipse = x;
            yEllipse = y;
        end
        LMplane.x(iRGC,:) = x;
        LMplane.y(iRGC,:) = y;
        if (noFit)
             LMplane.xFit(iRGC,:) = nan(1,100);
             LMplane.yFit(iRGC,:) = nan(1,100);
        else
            LMplane.xFit(iRGC,:) = xEllipse;
            LMplane.yFit(iRGC,:) = yEllipse;
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
                    
    for targetRGC = 1:size(LMplane.x,1)
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
    
    figure(1);
    clf;
    for targetSFindex = 1:size(LMtuning,3)
        subplot(3,5,targetSFindex);
        hold on
        for cellIndex = 1:size(LMtuning,2)
            x = squeeze(LMtuning(:,cellIndex, targetSFindex)) .* cosd(LMangles(:));
            y = squeeze(LMtuning(:,cellIndex, targetSFindex)) .* sind(LMangles(:));
            % make it symmetric
            x = cat(1, x, -x);
            y = cat(1, y, -y);
            [xEllipse,  yEllipse, semiAxes, rfCenter, noFit] = fitEllipseToContour(x,y);
            plot(x,y,'ko');
            if (noFit)
                plot(x,y,'-');
            else
                plot(xEllipse, yEllipse, 'k-');
            end
            set(gca, 'XLim', max(LMtuning(:))*[-1 1], 'YLim', max(LMtuning(:))*[-1 1]);
            axis 'square';
            pause
        end
        
        
    end
    
    
end
