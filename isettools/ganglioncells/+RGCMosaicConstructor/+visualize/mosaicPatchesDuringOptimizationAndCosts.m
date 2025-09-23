function mosaicPatchesDuringOptimizationAndCosts(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, minConeWeightVisualized, pStruct, intermediateMetaDataStructs, theVisualizedIntermediateConnectivityStages, visualizeSourceLatticeInsteadOfConnectedRFcenters,spatialChromaticUniformityTradeoff)

	patchesNum = numel(pStruct.xo);

    contourGenerationMethod = 'ellipseFitToPooledConePositions';
    maxNumberOfConesOutsideContour = pStruct.maxNumberOfConesOutsideContour;

    theFullSaturationCostColors = brewermap(3, 'Pastel1');

	if (~isempty(intermediateMetaDataStructs))
        if (isempty(theVisualizedIntermediateConnectivityStages))
            theVisualizedIntermediateConnectivityStage = 3:numel(intermediateMetaDataStructs);
        end
        innerLoopCount = numel(theVisualizedIntermediateConnectivityStage);
    else
        innerLoopCount = 1;
        plotTitle = '';
    end

    for iPatch = 1:patchesNum
        if (~isempty(intermediateMetaDataStructs))
            hFigPerformance = figure(2000+iPatch); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
            thePerformanceAxes = PublicationReadyPlotLib.generatePanelAxes(hFigPerformance,ff);
            hold(thePerformanceAxes{1,1}, 'on');

            hFigTransfers = figure(3000+iPatch); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theTransfersAxes = PublicationReadyPlotLib.generatePanelAxes(hFigTransfers,ff);
            hold(theTransfersAxes{1,1}, 'on');
        end % (~isempty(intermediateMetaDataStructs))

        allTransfers = [];
        allCosts = [];
        theCostNames = {};
        plotHandles = [];
        for iStage = 1:innerLoopCount  

            if (iStage == 1)
                theStageMarker = 'o';
            elseif (iStage == 2)
                theStageMarker = 'o';
            else
                theStageMarker = 'o'
            end
            costColorSaturation =  0.5*(innerLoopCount-iStage)/innerLoopCount;
            theCostColors = theFullSaturationCostColors*(1-costColorSaturation) + [0 0 0]*(costColorSaturation);

            % Reload the mosaic so we can crop it (to accelerate plotting)
            theImportedData = load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
            thePatchCenterDegs = [pStruct.xo(iPatch) pStruct.yo(iPatch)];

            if (~isempty(intermediateMetaDataStructs))
                dStruct = intermediateMetaDataStructs{theVisualizedIntermediateConnectivityStage(iStage)}
                theImportedData.theMRGCMosaic.debug_overwritePositionConnectivityData(dStruct);
            end

            % Crop it
            theCroppedMRGCMosaic = theImportedData.theMRGCMosaic.cropToSizeAtEccentricity(...
                    2*[pStruct.halfWidth pStruct.halfHeight], thePatchCenterDegs);
            clear 'theImportedData';

            fprintf('The cropped mosaic contains %d mRGCs\n', theCroppedMRGCMosaic.rgcsNum);

            % Compute spatial and chromatic uniformity
            [theSpatialCompactnessCosts, theSpectralUniformityCosts, ~, ...
            theInputConeNumerosityDifferentials, ...
            theCentroidOverlapCosts, theSpatialVarianceCosts] = theCroppedMRGCMosaic.rfCenterSpatioChromaticCosts();
        
            if (~isempty(intermediateMetaDataStructs))
                if (~isempty(dStruct.netTransfers))
                    transfersDuringPassesAtThisStage = dStruct.netTransfers;
                    costsDuringPassesAtThisStage = dStruct.netCostSequences;
                    if (isempty(allTransfers))
                        index0 = 0;
                    end
                    allTransfers = cat(1, allTransfers(:), transfersDuringPassesAtThisStage(:));
                    allCosts = cat(1, allCosts, costsDuringPassesAtThisStage);
                    if (contains(dStruct.phaseDescriptor, 'transfering'))
                        barFaceColor = [0.5 0.5 0.5];
                    else
                        barFaceColor = [0.8 0.8 0.8];
                    end
                    indices = index0 + (1:numel(transfersDuringPassesAtThisStage));

                    
                    % The spatial compactness cost
                    theSpatialCompactnessComponent1CostIndex = find(strcmp(dStruct.costComponentNames, 'numerosity differential cost'));
                    theSpatialCompactnessComponent2CostIndex = find(strcmp(dStruct.costComponentNames, 'spatial overlap cost'));
                    spatialCompactnessCosts = allCosts(indices, theSpatialCompactnessComponent1CostIndex) + ...
                                              allCosts(indices, theSpatialCompactnessComponent2CostIndex);
                    plotHandles(numel(plotHandles)+1)  = plot(thePerformanceAxes{1,1}, indices, spatialCompactnessCosts, ...
                        'k-', 'MarkerFaceColor', theCostColors(2,:), 'MarkerEdgeColor', [0.3 0.3 0.3],  ...
                        'Marker', theStageMarker, 'MarkerSize', ff.markerSize-2, 'LineWidth', ff.lineWidth); 
                    theCostNames{numel(theCostNames)+1} = 'spatial compactness (total)';

                    % The numerosity differential components
                    plotHandles(numel(plotHandles)+1)  = plot(thePerformanceAxes{1,1}, indices, allCosts(indices, theSpatialCompactnessComponent1CostIndex), ...
                        'k-', 'MarkerFaceColor', theCostColors(2,:), 'MarkerEdgeColor', [0.3 0.3 0.3],  ...
                        'Marker', 'v', 'MarkerSize', ff.markerSize-2, 'LineWidth', ff.lineWidth); 
                    theCostNames{numel(theCostNames)+1} = 'spatial compactness (numerosity differential)';

                    plotHandles(numel(plotHandles)+1)  = plot(thePerformanceAxes{1,1}, indices, allCosts(indices, theSpatialCompactnessComponent2CostIndex), ...
                        'k--', 'MarkerFaceColor', theCostColors(2,:), 'MarkerEdgeColor', [0.3 0.3 0.3],  ...
                        'Marker', '^', 'MarkerSize', ff.markerSize-2, 'LineWidth', ff.lineWidth); 
                    theCostNames{numel(theCostNames)+1} = 'spatial compactness (overlap)';


                    % The spectral purity cost
                    theCostIndex = find(strcmp(dStruct.costComponentNames, 'spectral uniformity cost'))
                    plotHandles(numel(plotHandles)+1)  = plot(thePerformanceAxes{1,1}, indices, allCosts(indices,theCostIndex), ...
                        'k-', 'MarkerFaceColor', theCostColors(1,:), 'MarkerEdgeColor', [0.3 0.3 0.3],  ...
                        'Marker', theStageMarker, 'MarkerSize', ff.markerSize-2, 'LineWidth', ff.lineWidth); 
                    theCostNames{numel(theCostNames)+1} = sprintf('%s', strrep(dStruct.costComponentNames{theCostIndex}, 'cost', ''));

                    % The total cost
                    theCostIndex = find(strcmp(dStruct.costComponentNames, 'total cost'))
                    plotHandles(numel(plotHandles)+1) = plot(thePerformanceAxes{1,1}, indices, allCosts(indices,theCostIndex), ...
                        'k-', 'MarkerFaceColor', theCostColors(3,:), 'MarkerEdgeColor', [0.3 0.3 0.3],  ...
                        'Marker', theStageMarker, 'MarkerSize', ff.markerSize-2, 'LineWidth', ff.lineWidth); 
                    theCostNames{numel(theCostNames)+1} = sprintf('%s', strrep(dStruct.costComponentNames{theCostIndex}, 'cost', ''));

                    XLims = [1 55];
                    YLims = [0 2];
                    set(thePerformanceAxes{1,1}, 'YLim', YLims, 'YTick', 0:0.2:10, ...
                        'XLim', XLims, 'XTick', 0:5:100);
                    legend(thePerformanceAxes{1,1}, plotHandles(1:5), theCostNames{1:5}, 'Location', 'NorthOutside', 'NumColumns', 2,  'Box', 'off');
                    xlabel(thePerformanceAxes{1,1}, 'pass no');
                    ylabel(thePerformanceAxes{1,1}, 'cost');

                    PublicationReadyPlotLib.offsetAxes(thePerformanceAxes{1,1}, ff, XLims, YLims, ...
                        'keepXaxisSymmetric', true, ...
                        'keepYaxisSymmetric', true);
                    

                    YtickIncrement = 50;
					

                    YLims = [0 (floor(max(allTransfers)/100)+1)*100];
                    YLims = [0 300];
                    bar(theTransfersAxes{1,1}, indices, allTransfers(indices), 1, 'FaceColor', barFaceColor);
                    set(theTransfersAxes{1,1}, 'YLim', YLims, ...
                    	       'XLim', XLims, 'XTick', 0:5:100);
                    set(theTransfersAxes{1,1}, 'YTick', YLims(1):YtickIncrement:YLims(2));
                    legend(theTransfersAxes{1,1}, {' cone transfers ', ' cone swaps '}, 'Location', 'NorthOutside', 'NumColumns', 2, 'Box', 'off');
                    xlabel(theTransfersAxes{1,1}, 'pass no');
                    ylabel(theTransfersAxes{1,1}, 'no. of reassignments');

                    PublicationReadyPlotLib.offsetAxes(theTransfersAxes{1,1}, ff, XLims, YLims, ...
                        'keepXaxisSymmetric', true, ...
                        'keepYaxisSymmetric', true);

                    index0 = indices(end);
                end

                plotTitle = sprintf('%s\n\\phi = %2.2f \\rightarrow (\\chi = %2.3f, \\lambda = %2.3f)', ...
                   dStruct.phaseDescriptor, spatialChromaticUniformityTradeoff, ...
                   mean(theSpatialCompactnessCosts), mean(theSpectralUniformityCosts));
            end % (~isempty(intermediateMetaDataStructs))

            

            domainVisualizationLimits = [...
                    pStruct.xo(iPatch) - 0.5*pStruct.visualizedWidth ...
                    pStruct.xo(iPatch) + 0.5*pStruct.visualizedWidth ...
                    pStruct.yo(iPatch) - 0.5*pStruct.visualizedHeight ...
                    pStruct.yo(iPatch) + 0.5*pStruct.visualizedHeight];
            domainVisualizationTicks = struct('x', -50:pStruct.tickIncrementDegs:50, 'y', -50:pStruct.tickIncrementDegs:50);
    
            thePatchSizeDegs = [0.1 0.1];
            theRGCindices = theCroppedMRGCMosaic.indicesOfRGCsWithinROI(thePatchCenterDegs, thePatchSizeDegs);
            theSelectRGCindex = theRGCindices(1);

            hFig = figure(1000+iPatch); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            
            if (visualizeSourceLatticeInsteadOfConnectedRFcenters)
                plotTitle = 'lattices';
            end

            theCroppedMRGCMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', theAxes{1,1}, ...
                    'domainVisualizationLimits', domainVisualizationLimits, ...
                    'domainVisualizationTicks', domainVisualizationTicks, ...
                    'minConeWeightVisualized', exp(-0.5), ...
                    'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
                    'identifiedConeApertureThetaSamples', 30, ...
                    'renderSourceLatticeInsteadOfConnectedRFcenters', visualizeSourceLatticeInsteadOfConnectedRFcenters, ...
                    'identifyInputCones', true, ...
                    'identifyPooledCones', true, ...
                    'inputConesAlpha', 0.9, ...
                    'pooledConesLineWidth', 1.5, ...
                    'centerSubregionContourSamples', 32, ...
                    'contourGenerationMethod', contourGenerationMethod, ...
                    'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
                    'plottedRFoutlineLineWidth', 1.0, ...
                    'plottedRFoutlineFaceAlpha', 0.5, ...
                    'plotTitle', plotTitle);
    
            % Finalize figure using the Publication-Ready format
            ff.box = 'on';
            PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
            
            % Export figure
            theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('_PatchXpos%2.3fDegs.pdf', pStruct.xo(iPatch)) ));
            NicePlot.exportFigToPDF(thePDFfileName, hFig,  300);

            if (~isempty(intermediateMetaDataStructs))
                thePDFfileNameForThisIntermediateStage = strrep(thePDFfileName, 'Degs.pdf', sprintf('Degs_%s.pdf', dStruct.phaseDescriptor));
                NicePlot.exportFigToPDF(thePDFfileNameForThisIntermediateStage, hFig,  300);
            end
        end % for iStage

        if (~isempty(intermediateMetaDataStructs))
        	ff.box = 'off';
            PublicationReadyPlotLib.applyFormat(thePerformanceAxes{1,1},ff);
            thePDFfileName = strrep(thePDFfileName, 'Degs', 'Degs_Performance');
            NicePlot.exportFigToPDF(thePDFfileName, hFigPerformance,  300);

            PublicationReadyPlotLib.applyFormat(theTransfersAxes{1,1},ff);
            thePDFfileName = strrep(thePDFfileName, 'Degs_Performance', 'Degs_Transfers');
            NicePlot.exportFigToPDF(thePDFfileName, hFigTransfers,  300);
        end

    end  % iPatch
end
