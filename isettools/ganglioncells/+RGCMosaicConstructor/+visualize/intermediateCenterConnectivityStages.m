function intermediateCenterConnectivityStages(theCenterConnectedMRGCMosaicFullFileName, ...
    theCenterConnectedMRGCMosaicFileName, intermediateConnectivityStageMetaDataFile, ...
    patchCoords, spatialChromaticUniformityTradeoff, varargin)
    

    % Parse input
    p = inputParser;
    p.addParameter('onlyShowRFoutlines', false, @islogical);
    p.addParameter('debugStage1', false, @islogical);
    p.parse(varargin{:});
    onlyShowRFoutlines = p.Results.onlyShowRFoutlines;
    debugStage1 = p.Results.debugStage1;


    % Load the intermediate meta data structs
    load(intermediateConnectivityStageMetaDataFile, 'intermediateMetaDataStructs');

    stagesNum = numel(intermediateMetaDataStructs);
    for iStage = 1:stagesNum
        visualizeConnectivityAtStage(theCenterConnectedMRGCMosaicFullFileName, ...
            theCenterConnectedMRGCMosaicFileName, intermediateMetaDataStructs{iStage}, ...
            iStage, patchCoords, spatialChromaticUniformityTradeoff, onlyShowRFoutlines,  debugStage1);
    end

end

function visualizeConnectivityAtStage(theCenterConnectedMRGCMosaicFullFileName, ...
    theCenterConnectedMRGCMosaicFileName, dStruct, iStage, patchCoords, ...
    spatialChromaticUniformityTradeoff, onlyShowRFoutlines, debugStage1)

    
    identifyInputCones = true; identifyPooledCones = true;
    contourGenerationMethod = 'ellipseFitToPooledConeApertureImage';

    if (iStage == 1)
        identifyPooledCones = false;
        contourGenerationMethod = 'ellipseFitBasedOnLocalSpacing';
    end

    plotTitle = dStruct.phaseDescriptor;
    if (iStage > 2)   
        plotTitleBase = sprintf('%s-%d, \\phi = %2.2f', plotTitle,iStage-2, spatialChromaticUniformityTradeoff);
    end

    % Visualize small patches of the mosaic along the horizontal meridian
    patchesNum = numel(patchCoords.xo);
    for iPatch = 1:patchesNum
                    
        % Reload the mosaic so we can crop it (to accelerate plotting)
        theImportedData = load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

        % Overwrite the connectivity data so we can visualize state at this stage
        theImportedData.theMRGCMosaic.debug_overwritePositionConnectivityData(dStruct);

        % Crop mosaic to accelerate plotting
        if (iStage == 1) && (debugStage1)
            theCroppedMRGCMosaic = ...
                theImportedData.theMRGCMosaic.cropToSizeAtEccentricity(10*[1 1], [patchCoords.xo(iPatch) patchCoords.yo(iPatch)]);
        else
            theCroppedMRGCMosaic = ...
                theImportedData.theMRGCMosaic.cropToSizeAtEccentricity(2*patchCoords.halfWidth*[1 1], [patchCoords.xo(iPatch) patchCoords.yo(iPatch)]);
        end


        if (iStage>2)
            % Compute spatial and chromatic uniformity
            [theSpatialUniformities, theChromaticUniformities] = theCroppedMRGCMosaic.rfCenterSpatioChromaticUniformity();
            plotTitle = sprintf('%s  \\rightarrow (\\chi = %2.2f, \\lambda = %2.2f)', plotTitleBase, mean(theSpatialUniformities), mean(theChromaticUniformities));
        end

        clear 'theImportedData';

        if (iStage == 1) && (debugStage1)
            domainVisualizationLimits = [...
                    patchCoords.xo(iPatch) - 5 ...
                    patchCoords.xo(iPatch) + 5 ...
                    patchCoords.yo(iPatch) - 5 ...
                    patchCoords.yo(iPatch) + 5];    
                    domainVisualizationTicks = struct('x', -50:1:50, 'y', -50:1:50);
        
        else
            domainVisualizationLimits = [...
                    patchCoords.xo(iPatch) - 0.5*patchCoords.visualizedWidth ...
                    patchCoords.xo(iPatch) + 0.5*patchCoords.visualizedWidth ...
                    patchCoords.yo(iPatch) - 0.5*patchCoords.visualizedHeight ...
                    patchCoords.yo(iPatch) + 0.5*patchCoords.visualizedHeight];   
            domainVisualizationTicks = struct('x', -50:patchCoords.tickIncrementDegs:50, 'y', -50:patchCoords.tickIncrementDegs:50);
        end

        

        hFig = renderIt(2000+iPatch, theCroppedMRGCMosaic, onlyShowRFoutlines, ...
            identifyInputCones, identifyPooledCones, contourGenerationMethod, plotTitle, ...
            domainVisualizationLimits, domainVisualizationTicks);

        % Export figure
        theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        if (onlyShowRFoutlines)
            thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('_xo%2.2f_atCenterConnectivityStage_%d_outlinesOnly.pdf',patchCoords.xo(iPatch), iStage) ));
        else
            thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('_xo%2.2f_atCenterConnectivityStage_%d.pdf',patchCoords.xo(iPatch), iStage) ));
        end

        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
     end % iPatch
end


function hFig = renderIt(figNo, theCroppedMRGCMosaic, onlyShowRFoutlines,  ...
    identifyInputCones, identifyPooledCones, contourGenerationMethod, plotTitle, ...
    domainVisualizationLimits, domainVisualizationTicks)
    hFig = figure(figNo); clf;
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    if (onlyShowRFoutlines)
        theCroppedMRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', theAxes{1,1}, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'centerSubregionContourSamples', 40, ...
            'contourGenerationMethod', contourGenerationMethod, ...
            'identifyInputCones', false, ...
            'identifyPooledCones', false, ...
            'plottedRFoutlineLineWidth', 1.0, ...
            'plotTitle', plotTitle);
    else
        theCroppedMRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', theAxes{1,1}, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'centerSubregionContourSamples', 40, ...
            'contourGenerationMethod', contourGenerationMethod, ...
            'identifyInputCones', identifyInputCones, ...
            'identifyPooledCones', identifyPooledCones, ...
            'pooledConesLineWidth', 1.5, ...
            'plottedRFoutlineLineWidth', 1.0, ...
            'plotTitle', plotTitle);
    end


    % Finalize figure using the Publication-Ready format
    ff.box = 'on';
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

end
