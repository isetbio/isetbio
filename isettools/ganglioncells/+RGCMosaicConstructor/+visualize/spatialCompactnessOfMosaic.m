function spatialCompactnessOfMosaic(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, minConeWeightVisualized, pStruct)

	contourGenerationMethod = 'ellipseFitToPooledConePositions';
	maxNumberOfConesOutsideContour = pStruct.maxNumberOfConesOutsideContour;

	nZones = pStruct.maxConesNum - pStruct.minConesNum+1;
    coneNumerosityLimits = zeros(nZones,2);
    coneNumerosityZones = pStruct.minConesNum + (0:(nZones-1));
    coneNumerosityTicks = coneNumerosityZones;


    if (nZones <= 5)
    	skippedTicks = 0;
    elseif (nZones <= 10)
    	skippedTicks = 1;
    elseif (nZones <= 20)
    	skippedTicks = 2;
    elseif (nZones <= 40)
    	skippedTicks = 5;
    else
    	skippedTicks = 10;  %max([1 ceil(nZones/8)])
    end


    for iZone = 1:nZones
        coneNumerosityLimits(iZone,:) = coneNumerosityZones(iZone)*[1 1+eps];
        if (skippedTicks>1)
	        if (mod(iZone,skippedTicks) == 1)
	        	coneNumerosityTickLabels{iZone} = sprintf('%d', coneNumerosityTicks(iZone));
	        else
	        	coneNumerosityTickLabels{iZone} = '';
	        end
	    else
	    	coneNumerosityTickLabels{iZone} = sprintf('%d', coneNumerosityTicks(iZone));
	    end
    end
    

    nSpatialVarianceZones = 33;
    spatialVarianceZones = linspace(0, 3, nSpatialVarianceZones);
    spatialVarianceFaceColors = brewermap(nSpatialVarianceZones, 'blues');
    spatialVarianceTicks = spatialVarianceZones;


    for iZone = 1:nSpatialVarianceZones
    	if (mod(iZone,8) == 1)
    	 	spatialVarianceTickLabels{iZone} = sprintf('%2.2f', spatialVarianceTicks(iZone));
    	 else
    	 	spatialVarianceTickLabels{iZone} = '';
    	 end
    end

    nSpatialSeparationZones = 33;
    spatialSeparationZones = linspace(0,3, nSpatialSeparationZones);
    spatialSeparationFaceColors = brewermap(nSpatialVarianceZones, 'blues');
    spatialSeparationTicks = spatialSeparationZones;
    for iZone = 1:nSpatialSeparationZones
    	if (mod(iZone,8) == 1)
    	 	spatialSeparationTickLabels{iZone} = sprintf('%2.2f', spatialSeparationTicks(iZone));
    	 else
    	 	spatialSeparationTickLabels{iZone} = '';
    	 end
    end


    plottedRFoutlineFaceColors = brewermap(nZones, 'spectral');

	for iXo = 1:numel(pStruct.rapheXo)
		bandWidth = pStruct.rapheBandWidth;
       	bandHeight = pStruct.rapheBandHeight;
       	bandXo = pStruct.rapheXo(iXo);
        bandYo = pStruct.rapheYo(iXo);

		% Reload the mosaic so we can crop it (to accelerate plotting)
        theImportedData = load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

        % Crop it
        theCroppedMRGCMosaic = theImportedData.theMRGCMosaic.cropToSizeAtEccentricity([bandWidth bandHeight], [bandXo bandYo]);
        clear 'theImportedData';

        % Visualized domain
        xyMin = min(theCroppedMRGCMosaic.rgcRFpositionsDegs,[],1);
        xyMax = max(theCroppedMRGCMosaic.rgcRFpositionsDegs,[],1);
        domainVisualizationLimits = [bandXo-bandWidth/2 bandXo+bandWidth/2 bandYo-bandHeight/2 bandYo+bandHeight/2];
        domainVisualizationTicks = struct('x', -50:1:50, 'y', -50:1:50);
 	

        hFig = figure(1000); clf;
        ff = PublicationReadyPlotLib.figureComponents(pStruct.panelFormat);
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

        % Compute thePatchSpectralUniformity, thePatchCentroidOverlap, thePatchCentroidSigma
        [theSpatialCompactnessCosts, theSpectralUniformityCosts, ~, ~, ...
         theCentroidOverlapCosts, theSpatialVarianceCosts] = theCroppedMRGCMosaic.rfCenterSpatioChromaticCosts();

        theRGCConeMixtures = 0.5 + (0.5-theSpectralUniformityCosts/2);
		theRGCCentroidSeparations = 1./theCentroidOverlapCosts;
		theRGCCentroidSigmas = theSpatialVarianceCosts;
		
		% Compute histogram of spatial variances of centroids
		spatialVarianceCounts = histcounts(theRGCCentroidSigmas, spatialVarianceZones);
		spatialVarianceCounts = spatialVarianceCounts / theCroppedMRGCMosaic.rgcsNum;

		% Compute histogram of spatial separation of centroids
		spatialSeparationCounts = histcounts(theRGCCentroidSeparations, spatialSeparationZones);
		spatialSeparationCounts = spatialSeparationCounts / theCroppedMRGCMosaic.rgcsNum;


		plotTitle = '';
		coneNumerosityHistogram = zeros(1, numel(coneNumerosityZones));
		
		for iZone = numel(coneNumerosityZones):-1:1
			targetConeNumerosity = coneNumerosityZones(iZone);
			targetRGCindices = theCroppedMRGCMosaic.indicesOfRGCsWithTargetCenterConesNumInRange(...
                        coneNumerosityLimits(iZone,:), ...
                        'minConeWeightIncluded', minConeWeightVisualized);

            coneNumerosityHistogram(iZone) = numel(targetRGCindices)/theCroppedMRGCMosaic.rgcsNum;
            if (numel(targetRGCindices)>0)
            	theCroppedMRGCMosaic.visualize(...
	                'figureHandle', hFig, ...
	                'axesHandle', theAxes{1,1}, ...
	                'domainVisualizationLimits', domainVisualizationLimits, ...
	                'domainVisualizationTicks', domainVisualizationTicks, ...
	                'visualizedRGCindices', targetRGCindices, ...
	                'plottedRFoutlineFaceColor', plottedRFoutlineFaceColors(iZone,:), ...
	                'plottedRFoutlineFaceAlpha', 0.75, ...
	                'minConeWeightVisualized', minConeWeightVisualized, ...
	                'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
	                'identifyPooledCones', ~true, ...
	                'identifyInputCones', ~true, ...
	                'contourGenerationMethod', contourGenerationMethod, ...
	                'clearAxesBeforeDrawing', false, ...
	                'centerSubregionContourSamples', 20, ...
	                'plottedRFoutlineLineWidth', 1.0, ...
	                'plotTitle', plotTitle);
            end % if (numel(targetRGCindices)>0)

            % Set the colormap
            colormap(theAxes{1,1},plottedRFoutlineFaceColors);
            set(theAxes{1,1}, 'CLim', [0 size(plottedRFoutlineFaceColors,1)]);

            % Finalize figure using the Publication-Ready format
        	ff.box = 'on';
        	PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
		end % iZone

        % Add the colorbar
        c = colorbar(theAxes{1,1}, ...
        			'Location', 'NorthOutside', ...
                    'Ticks', coneNumerosityTicks-0.5, ...
                    'TickLabels', coneNumerosityTickLabels);
        c.TickLength = 0.0;
        c.Label.String = '# of cones';

        % Finalize figure using the Publication-Ready format
        ff.box = 'on';
        PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
                 
        if (pStruct.includeHistogramInsets) 
	        % Add the cone numerosity histogram
	        histogramAxes = axes(hFig, 'Position',[0.05 0.10 0.15 0.2]);
	        b = bar(histogramAxes, coneNumerosityTicks, coneNumerosityHistogram, 1);
	        b.FaceColor = 'flat';
	        for iBar = 1:size(plottedRFoutlineFaceColors,1)
				b.CData(iBar,:) = plottedRFoutlineFaceColors(iBar,:);
			end
			XLims = [min(coneNumerosityTicks)-0.5 max(coneNumerosityTicks)+0.5];
			YLims = [0 1];
			set(histogramAxes, 'YLim', YLims, 'YTick', 0:0.2:1, 'XLim', XLims, ...
				'XTick', coneNumerosityTicks, 'XTickLabels', coneNumerosityTickLabels);

			grid(histogramAxes, 'on');
			box(histogramAxes, 'on');
			PublicationReadyPlotLib.offsetAxes(histogramAxes, ff, XLims, YLims, ...
				'keepXaxisSymmetric', true);
			PublicationReadyPlotLib.applyFormat(histogramAxes,ff);
		  	xlabel(histogramAxes, '# of input cones')
		  	ylabel(histogramAxes, 'fraction of total RGCs')


		  	% Add the spatial variance histogram
		  	histogramAxes = axes(hFig, 'Position',[0.05 0.42 0.15 0.2]);
		  	bb = bar(histogramAxes, spatialVarianceZones(1:end-1), spatialVarianceCounts,1);
		  	bb.FaceColor = 'flat';
	        for iBar = 1:(size(spatialVarianceFaceColors,1)-1)
				bb.CData(iBar,:) = spatialVarianceFaceColors(iBar,:);
			end
			XLims = [spatialVarianceZones(1) spatialVarianceZones(end)];
			YLims = [0 1];
			set(histogramAxes, 'YLim', YLims, 'YTick', 0:0.2:1, 'XLim', XLims, ...
				'XTick', spatialVarianceTicks, 'XTickLabels', spatialVarianceTickLabels);
			grid(histogramAxes, 'on');
			box(histogramAxes, 'on');
			PublicationReadyPlotLib.offsetAxes(histogramAxes, ff, XLims, YLims, ...
				'keepXaxisSymmetric', true);
			PublicationReadyPlotLib.applyFormat(histogramAxes,ff);
		  	xlabel(histogramAxes, 'spatial variance of input cones')
		  	ylabel(histogramAxes, 'fraction of total RGCs')

		  	% Add the spatial separation histogram
		  	histogramAxes = axes(hFig, 'Position',[0.05 0.72 0.15 0.2]);
		  	bb = bar(histogramAxes, spatialSeparationZones(1:end-1), spatialSeparationCounts,1);
		  	bb.FaceColor = 'flat';
	        for iBar = 1:(size(spatialSeparationFaceColors,1)-1)
				bb.CData(iBar,:) = spatialSeparationFaceColors(iBar,:);
			end
			XLims = [spatialSeparationZones(1) spatialSeparationZones(end)];
			YLims = [0 1];
			set(histogramAxes, 'YLim', YLims, 'YTick', 0:0.2:1, 'XLim', XLims, ...
				'XTick', spatialSeparationTicks, 'XTickLabels', spatialSeparationTickLabels);
			grid(histogramAxes, 'on');
			box(histogramAxes, 'on');
			PublicationReadyPlotLib.offsetAxes(histogramAxes, ff, XLims, YLims, ...
				'keepXaxisSymmetric', true);
			PublicationReadyPlotLib.applyFormat(histogramAxes,ff);
		  	xlabel(histogramAxes, 'separation')
		  	ylabel(histogramAxes, 'fraction of total RGCs');
		end % if (pStruct.includeHistogramInsets) 



        % Export figure
        theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('_ConeNumerosityBand_%2.1fDegs_%2.1fDegs.pdf',bandXo, bandYo) ));
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
   	end % iXo
end
