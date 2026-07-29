function spectralUniformityOfMosaic(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, minConeWeightVisualized, pStruct)

	contourGenerationMethod = 'ellipseFitToPooledConePositions';

	maxNumberOfConesOutsideContour = 2;

	coneMixtureZones = [0.5 2/3 3/4 5/6 1];
	coneMixtureLimits(1,:) = [0.5 0.6];
	coneMixtureLimits(2,:) = [0.6 0.7];
	coneMixtureLimits(3,:) = [0.7 0.8];
	coneMixtureLimits(4,:) = [0.8 0.9];
	coneMixtureLimits(5,:) = [0.9 1.0+eps];

	coneMixtureTickLabels = {'1/2', '2/3', '3/4', '5/6', '1/1'};
	coneMixtureTicks = (1:numel(coneMixtureZones)) - 0.5;

	plottedRFoutlineFaceColors = brewermap(numel(coneMixtureZones), '*spectral');

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
		theRGCCentroidOverlaps = 1./theCentroidOverlapCosts;
		theRGCCentroidSigmas = theSpatialVarianceCosts;

		plotTitle = '';
		coneMixtureHistogram = zeros(1, numel(coneMixtureZones));
		for iZone = numel(coneMixtureZones):-1:1
			targetConeMixture = coneMixtureZones(iZone);
            targetRGCindices = find(theRGCConeMixtures>=coneMixtureLimits(iZone,1) & theRGCConeMixtures<coneMixtureLimits(iZone,2));
            [iZone numel(targetRGCindices)]
            coneMixtureHistogram(iZone) = numel(targetRGCindices)/theCroppedMRGCMosaic.rgcsNum;
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
                    'Ticks', coneMixtureTicks, ...
                    'TickLabels', coneMixtureTickLabels);
        c.TickLength = 0.0;
        c.Label.String = 'cone mixture';

        % Finalize figure using the Publication-Ready format
        ff.box = 'on';
        PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
                   
        % Add the histogram
        if (contains(pStruct.panelFormat, 'standard'))
        	histogramAxesPosition = [0.02 0.02 0.2 0.25];
        else
        	histogramAxesPosition = [0.05 0.10 0.15 0.2];
        end

        histogramAxes = axes(hFig, 'Position', histogramAxesPosition);
        b = bar(histogramAxes, coneMixtureTicks, coneMixtureHistogram, 1);
        b.FaceColor = 'flat';
        for iBar = 1:size(plottedRFoutlineFaceColors,1)
			b.CData(iBar,:) = plottedRFoutlineFaceColors(iBar,:);
		end
		XLims = [0.0 numel(coneMixtureTicks)];
		YLims = [0 1];

		if (contains(pStruct.panelFormat, 'standard'))
			set(histogramAxes, 'YLim', YLims, 'YTick', 0:0.2:1, 'XLim', XLims, ...
				'XTick', coneMixtureTicks, 'XTickLabels', {}, 'YTickLabels', {});
		else
			set(histogramAxes, 'YLim', YLims, 'YTick', 0:0.2:1, 'XLim', XLims, ...
				'XTick', coneMixtureTicks, 'XTickLabels', coneMixtureTickLabels);
		end
		grid(histogramAxes, 'on');
		box(histogramAxes, 'on');
		PublicationReadyPlotLib.offsetAxes(histogramAxes, ff, XLims, YLims, ...
			'keepXaxisSymmetric', true);
		PublicationReadyPlotLib.applyFormat(histogramAxes,ff);
		if (~contains(pStruct.panelFormat, 'standard'))
	  		xlabel(histogramAxes, 'cone mixture')
	  		ylabel(histogramAxes, 'fraction of total RGCs')
	  	end

        % Export figure
        theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('_ChromaUniformityBand_%2.1fDegs_%2.1fDegs.pdf',bandXo, bandYo) ));
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
   	end % iXo
end
