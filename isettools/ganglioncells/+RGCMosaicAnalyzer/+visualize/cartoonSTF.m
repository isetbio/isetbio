function cartoonSTF()
%
%   RGCMosaicAnalyzer.visualize.cartoonSTF()
%
	% RGCMosaicAnalyzer.visualize.cartoonSTF
	

	hFig = figure(12); clf;
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	ax = theAxes{1,1};
	
	dataPointsSF = logspace(log10(0.1), log10(60), 16);
	%dataPointsSF = 3;
	showDataPoints = true;
	showCompositeSTF = true;
	showComponentSTFs = ~true;
	surroundStrengthModifier = 0.0;

	thePDFfileName = 'STF.pdf';
	generateSTF(hFig, ax, ff, dataPointsSF, showDataPoints, showCompositeSTF, showComponentSTFs, surroundStrengthModifier, thePDFfileName);


	hFig = figure(11); clf;
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
	ax = theAxes{1,1};

	STFamp = 0.4;
	frequency = 3;
	thePDFfileName = 'sinusoid.pdf';
	generateSinusoid(hFig, ax, ff, STFamp, frequency, thePDFfileName);



	hFig = figure(10); clf;
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};


    contrast = 0.9
    chromaticity = 'MconeIsolating';
    chromaticity = 'Achromatic';
    spatialFrequencyCPD = [4 1 8 2];

    theMovieFileName = sprintf('%sGrating', chromaticity);
	generateDriftinGrating(ax, spatialFrequencyCPD, contrast, chromaticity, theMovieFileName);
end

function generateDriftinGrating(ax,  spatialFrequencyCPD, contrast, chromaticity, theMovieFileName);

    XLims = [-2 2];
    YLims = [-2 2];
    
    x = linspace(-2.05,2.05,256);
    y = x;
    [X,Y] = meshgrid(x,y);
    sigmaX = 1.5;
    gaussianEnvelope = exp(-0.5*(X/sigmaX).^2) .* exp(-0.5*(Y/sigmaX).^2);

	[coneContrasts, totalContrast] = ...
    		visualStimulusGenerator.coneContrastsFromChromaticity(chromaticity);

    viewingDistanceMeters = 4;
    resolutionDegs = 0.01;
    wavelengthSupport = 400:10:700;
  	thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
            wavelengthSupport, resolutionDegs, ...
            viewingDistanceMeters);


    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = 'demos';

    theVideoFileName = sprintf('%s/%s/%s',theRawFiguresDir, pdfExportSubDir,theMovieFileName);
	videoOBJ = VideoWriter(theVideoFileName, 'MPEG-4'); % H264 format
	videoOBJ.FrameRate = 30;
	videoOBJ.Quality = 100;
	videoOBJ.open();

	for iSF = 1:numel(spatialFrequencyCPD)
    	stimParams = struct(...
		    'backgroundChromaticity', [0.301 0.301], ...
		    'backgroundLuminanceCdM2', 50, ...
		    'contrast', totalContrast, ...
		    'coneContrasts', coneContrasts, ...
		    'sizeDegs', 1, ...
		    'positionDegs', [0 0], ...
		    'resolutionDegs', 0.002, ...
		    'spatialPhaseIncrementDegs', 10, ...
		    'orientationDegs', 0, ...
		    'spatialFrequencyCPD', spatialFrequencyCPD(iSF), ...
		    'temporalFrequencyHz', 3.0, ...
		    'durationSeconds', 0.5, ...
		    'temporalEnvelopeTau', 0.15, ...
		    'coneMosaicModulationBasedResponse',  true ...
		    );

    
		[theDriftingGratingSpatialModulationPatterns, spatialSupportDegs, spatialPhasesDegs, ...
       		temporalSupportSeconds, temporalRamp] = visualStimulusGenerator.driftingGratingModulationPatterns(stimParams);

       	figure(99);
       	plot(temporalSupportSeconds, temporalRamp, 'k-')
       	drawnow;

    	[theDriftingGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
              thePresentationDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
              'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
              'validateScenes', ~true);

		for iPhase = 1:numel(theDriftingGratingFrameScenes)
	    	cla(ax);
	    	theSceneRGBimage = sceneGet(theDriftingGratingFrameScenes{iPhase}, 'rgbimage');

		    image(ax, spatialSupportDegs, spatialSupportDegs, theSceneRGBimage);
		    axis(ax, 'image');

		    set(ax, 'YLim', [min(spatialSupportDegs) max(spatialSupportDegs)], 'XLim',  [min(spatialSupportDegs) max(spatialSupportDegs)], ...
		    	'XTickLabel', {}, 'YTickLabel', {});
		    
		    grid(ax, 'off');
		    set(ax, 'XColor', 'none', 'YColor', 'none');
		    drawnow;
		    videoOBJ.writeVideo(getframe(ax));
		end
	end % iSF

	videoOBJ.close();

end



function generateSinusoid(hFig, ax, ff, STFamp, frequency, thePDFfileName)
	
    baseline = 0.0;
    faceColor = [1 0.5 0.5]; 
    edgeColor = [1 0 0];
    faceAlpha = [0.3]; 
    lineWidth = 1.5;
    x = linspace(-2.05,2.05,256);
    RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax, ...
    	x, STFamp * cosd(360*frequency*x), baseline, faceColor, edgeColor, faceAlpha, lineWidth);
    XLims = [-2 2];
    YLims = [-1 1];
    set(ax, 'CLim', [-1 1], 'YLim', YLims, 'XLim', XLims, 'XTick', -2:0.5:2, 'YTick', -1:0.5:1, 'XTickLabel', {}, 'YTickLabel', {});
	xlabel(ax, 'time (s)');
	ylabel(ax, 'response');

	% Finalize figure using the Publication-Ready format
	PublicationReadyPlotLib.applyFormat(ax,ff);
	%PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	if (~isempty(thePDFfileName))
		theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
		pdfExportSubDir = 'demos';
    	thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFfileName);
    	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end

function generateSTF(hFig, ax, ff, dataPointsSF, showDataPoints, showCompositeSTF, showComponentSTFs, surroundStrengthModifier, thePDFfileName)

    %                          Kc           intStoCsens             RsToRc            RcDegs    
    DoGparams.finalValues = [1.2313e+03     0.3207                  5.6952            0.0163];
    DoGparams.finalValues = [554.8430       0.9587                  8.0829            0.0132];
    DoGparams.finalValues = [701.0495       0.5770                  7.5501            0.0224];
    DoGparams.names         = {'Kc',        'intStoCsens',         'RsToRc',         'RcDegs'};
    
    DoGparams.finalValues(2) = surroundStrengthModifier*DoGparams.finalValues(2);
    % The DoG model in the frequency domain
    DoGSTF = @(params,spatialFrequency)(...
                    abs(params(1) * ( pi * params(4)^2 * exp(-(pi*params(4)*spatialFrequency).^2) ) - ...
                    params(1)*params(2)/(params(3))^2 * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*spatialFrequency).^2) )));
   
    DoGSTFcenter = @(params,spatialFrequency)(...
                    abs(params(1) * ( pi * params(4)^2 * exp(-(pi*params(4)*spatialFrequency).^2) )));

    DoGSTFsurround = @(params,spatialFrequency)(...
                    params(1) * params(2)/(params(3))^2 * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*spatialFrequency).^2) ));
   
    spatialFrequency = logspace(log10(0.1), log10(100), 64);

	%theVisualSTF = DoGSTF(DoGparams.finalValues,spatialFrequency);    
    theVisualSTFcenter = DoGSTFcenter(DoGparams.finalValues,spatialFrequency);
    theVisualSTFsurround = DoGSTFsurround(DoGparams.finalValues,spatialFrequency);
    theVisualSTF = theVisualSTFcenter - theVisualSTFsurround;

    if (showCompositeSTF)
	    plot(ax, spatialFrequency, theVisualSTF, 'k-', ...
	    	'Color', [1 0.5 0.5], ...
	    	'MarkerFaceColor', [1 0.5 0.5],...
	    	'LineWidth', 5.0);
	end
    hold(ax, 'on');
	if (showComponentSTFs)
	    plot(ax, spatialFrequency, theVisualSTFcenter, 'k-', ...
	    	'Color', [0.6 0.6 0.6], ...
	    	'LineWidth', 4.0);
	    plot(ax, spatialFrequency, theVisualSTFcenter, 'k--', ...
	    	'Color', [0.9 0.9 0.9], ...
	    	'LineWidth', 4.0);

	    plot(ax, spatialFrequency, theVisualSTFsurround, 'k-', ...
	    	'Color', [0.6 0.6 0.6], ...
	    	'LineWidth', 4.0);
	    plot(ax, spatialFrequency, theVisualSTFsurround, 'k--', ...
	    	'Color', [0.3 0.3 0.3], ...
	    	'LineWidth', 4.0);
	end

    if (showDataPoints)
    	for iPoint = 1:numel(dataPointsSF)
		    [~,idx] = min(abs(spatialFrequency-dataPointsSF(iPoint)));
		    plot(ax, spatialFrequency(idx), theVisualSTF(idx), 'ro-', ...
		    	'MarkerFaceColor', [1 0.5 0.5],...
		    	'MarkerSize', ff.markerSize+8, 'LineWidth', 2.0);
		end
	end

    XLims = [0.1 100];
    YLims = [0 max(theVisualSTF(:))*1.1];
    set(ax, 'XScale', 'log');
    set(ax, 'YLim', YLims, 'XLim', XLims, 'XTick', [0.1 0.3 1 3 10 30 100], ...
    	'YTick', linspace(YLims(1), YLims(2), 5), 'YTickLabel', {});
    
	xlabel(ax, 'spatial frequency (c/deg)');
	ylabel(ax, 'modulation amplitude');
	% Finalize figure using the Publication-Ready format
	PublicationReadyPlotLib.applyFormat(ax,ff);
	%PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

	if (~isempty(thePDFfileName))
		theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
		pdfExportSubDir = 'demos';
    	thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFfileName);
    	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end

end
