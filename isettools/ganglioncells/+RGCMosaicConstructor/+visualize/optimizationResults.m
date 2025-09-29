function optimizationResults(optimizationResultsStruct, ...
	theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, ...
	pdfExportSubDir, targetRFcenterConesNum, targetRFcenterDominantConeType,...
    summaryFigureInsteadOfSeparateFigures, theSummaryPDFfileName, varargin)

    p = inputParser;
    p.addParameter('fixedSpatialSupportTickSeparationArcMin', @(x)(isempty(x)||isscalar(x)));
    p.addParameter('fixedScalaBarDegs', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConeApertureImage', ...
        @(x)(ismember(x, mRGCMosaic.validRFsubregionContourGenerationMethods)));
    p.addParameter('maxNumberOfConesOutsideContour', 1, @isscalar);
    p.parse(varargin{:});

    fixedSpatialSupportTickSeparationArcMin = p.Results.fixedSpatialSupportTickSeparationArcMin;
    fixedScalaBarDegs = p.Results.fixedScalaBarDegs;
    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;
    contourGenerationMethod = p.Results.contourGenerationMethod;

    axSTF = [];
    axDoGFitParams = [];
    axSurroundFitParams = [];
    axCorrespondenceToPackerDaceyH1 = [];
    axAchievedPerformance = [];
    axCenterSurroundCompositeConePoolingMap = [];
    axSurroundConePoolingMap = [];
    axSurroundOnlyConePoolingMap = [];
    axHorizontalLineWeightingFunctions = [];
    axVerticalLineWeightingFunctions = [];

    if (summaryFigureInsteadOfSeparateFigures)
        hFigSummary = figure(1000); clf;
        set(hFigSummary, 'Name', 'Summary', 'Position', [10 10 2000 1000], 'Color', [1 1 1]);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 4, ...
           'heightMargin',  0.16, ...
           'widthMargin',    0.06, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.10, ...
           'topMargin',      0.02);

        axSTF = subplot('Position', subplotPosVectors(1,1).v);
        axDoGFitParams = subplot('Position', subplotPosVectors(1,2).v);
        axSurroundFitParams = subplot('Position', subplotPosVectors(1,3).v);
        axCorrespondenceToPackerDaceyH1 = subplot('Position', subplotPosVectors(1,4).v);
        axAchievedPerformance = subplot('Position', subplotPosVectors(2,1).v);
        axCenterSurroundCompositeConePoolingMap = subplot('Position', subplotPosVectors(2,2).v);
        axHorizontalLineWeightingFunctions = subplot('Position', subplotPosVectors(2,3).v);
        axVerticalLineWeightingFunctions = subplot('Position', subplotPosVectors(2,4).v);
    end


    % The STF
	figNo = 1; figPos = [50 700];

	renderSTF(pdfExportSubDir, figNo, figPos, ...
		optimizationResultsStruct.theFinalSTFdata.spatialFrequencySupport, ...
		optimizationResultsStruct.theFinalSTFdata.visualSTF, ...
		optimizationResultsStruct.theFinalSTFdata.fittedDoGModelToVisualSTF.sfHiRes, ...
		optimizationResultsStruct.theFinalSTFdata.fittedDoGModelToVisualSTF.compositeSTFHiRes, ...
		optimizationResultsStruct.theFinalSTFdata.fittedDoGModelToVisualSTF.centerSTFHiRes, ...
		optimizationResultsStruct.theFinalSTFdata.fittedDoGModelToVisualSTF.surroundSTFHiRes, ...
        'axesToRenderIn', axSTF);

	% The DoG fit to the cell's optimized STF (matching Rs/Rc, intS/C to C&K '95')
	figNo = 2; figPos = [700 700];
	renderModelParamsAndRanges(pdfExportSubDir, figNo, figPos, ...
		optimizationResultsStruct.theFinalSTFdata.fittedDoGModelParams, 'DoG fit to STF', ...
        'axesToRenderIn', axDoGFitParams);

	% The cell's optimized surround pooling model params (which leads to the above optimized STF, i.e., matching Rs/Rc, intS/C to C&K '95')
	figNo = 3; figPos = [1400 700];
    if (strcmp(optimizationResultsStruct.modelConstants.poolingModel.name, 'PackerDacey2002H1FixedCellIndex'))
        H1cellString = sprintf('%d', optimizationResultsStruct.modelConstants.poolingModel.fixedH1CellIndex);
        modelDescriptorString = sprintf('%s (H1-%s)', optimizationResultsStruct.modelConstants.poolingModel.name, H1cellString);
      else
        H1cellString = 'free';
        modelDescriptorString = sprintf('%s', optimizationResultsStruct.modelConstants.poolingModel.name);
      end

	renderModelParamsAndRanges(pdfExportSubDir, figNo, figPos, ...
		optimizationResultsStruct.modelVariables, modelDescriptorString, ...
        'axesToRenderIn', axSurroundFitParams);

	% The cell's fitted H1 params and their correspondence to Packer & Dacey data
	figNo = 4; figPos = [50 120];
	renderCorrespondenceBetweenFittedH1paramsAndPackerDaceyData(pdfExportSubDir, figNo, figPos, ...
 		optimizationResultsStruct.modelVariables, optimizationResultsStruct.modelConstants.poolingModel, ...
        'axesToRenderIn', axCorrespondenceToPackerDaceyH1);


    % The achieved ratios
    figNo = 5; figPos = [700 120];

    % Visualize the fit performance
    renderAchievedOptimizationPerformance(pdfExportSubDir, figNo, figPos, ...
        targetVisualSTFparams, optimizationResultsStruct.theFinalSTFdata.fittedDoGModelParams, ...
        optimizationResultsStruct.theFinalSTFdata.residual, ...
        optimizationResultsStruct.theFinalRMSE, ...
        H1cellString, ...
        targetRFcenterConesNum, targetRFcenterDominantConeType, ...
        'axesToRenderIn', axAchievedPerformance);

    % The cell's optimized pooling weights map.
    % Center pooling weights at exp(-0.5), surround pooling weights at 5%
    figNo = 8; figPos = [175 460]; 

    % Criterion used by Field (2010) for center cones:
    % (from supplemental methods section titled: "Cones providing input to receptive field center and surround"
    % "Cones in the RF center were defined as those with a weight greater than or equal to 10% of the weight of cone 
    % with the largest weight (which is by definition in the center)."
    % minCenterConeWeight = exp(-0.5);
    minCenterConeWeight = 0.1;

    % Criterion used by Field (2010) for surround cones:
    % (from supplemental methods section titled: "Cones providing input to receptive field center and surround"
    % "Surround cones were defined as those cones with a light response polarity opposite to the center, 
    % weights >0.5% that of the peak cone but of opposite sign, and within 8 SDs of a circular Gaussian fit to the RF center
    % Note however, that their measured center weights are composite, center-surround as you cant
    % separate the surround response from the center response electrophysiologically, so we should do the same here)
    minSurroundConeWeight = 0.5/100;

    % The spatial support
    spatialSupportCenterDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);
    ecc = sqrt(sum(spatialSupportCenterDegs.^2,2));

    if (ecc > 20)
        scaleBarDegs = 0.1;
        spatialSupportTickSeparationArcMin = 33.0;
    elseif (ecc > 10)
        scaleBarDegs = 0.08;
        spatialSupportTickSeparationArcMin = 16.0;
    elseif (ecc > 5)
        scaleBarDegs = 0.07;
        spatialSupportTickSeparationArcMin = 12.0;
    elseif (ecc > 2)
        scaleBarDegs = 0.06;
        spatialSupportTickSeparationArcMin = 8.0;
    else
        scaleBarDegs = 0.05;
        spatialSupportTickSeparationArcMin = 4.0;
    end

    if (~isempty(fixedSpatialSupportTickSeparationArcMin))
        spatialSupportTickSeparationArcMin = fixedSpatialSupportTickSeparationArcMin;
    end

    if (~isempty(fixedScalaBarDegs))
        scaleBarDegs = fixedScalaBarDegs;
    end


    if (~summaryFigureInsteadOfSeparateFigures)
        % surround alone (together with center ellipse)
        minSurroundConeWeightRelativity = 'center';
        visualizedConeWeights = 'surround-alone';
        plotTitle = sprintf('visualized cone weights: ''%s''\n Wc > %1.3f, Ws > %1.3f (x peak %s cone)',...
            visualizedConeWeights, minCenterConeWeight, minSurroundConeWeight, minSurroundConeWeightRelativity);

        RGCMosaicConstructor.visualize.conePoolingWeightsMap(...
            pdfExportSubDir, figNo, figPos, theMRGCMosaic, theTargetRGCindex, ...
            spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
            optimizationResultsStruct.theFinalPooledConeIndicesAndWeights, ...
            minCenterConeWeight, ...
            minSurroundConeWeight, ...
            minSurroundConeWeightRelativity, ...
            visualizedConeWeights, ...
            'contourGenerationMethod', contourGenerationMethod, ...
            'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
            'scaleBarDegs', scaleBarDegs, ...
            'plotTitle', plotTitle, ...
            'axesToRenderIn', axSurroundOnlyConePoolingMap);
    end

    % The cell's optimized pooling weights map. Just the RF center now
    figNo = 9; figPos = [175 560]; 

    % center-surround
    minSurroundConeWeightRelativity = 'center';
    visualizedConeWeights = 'center-surround';
    visualizedConeWeights = 'surround-alone';
    plotTitle = sprintf('visualized cone weights: ''%s''\n Wc > %1.3f, Ws > %1.3f (x peak %s cone)',...
        visualizedConeWeights, minCenterConeWeight, minSurroundConeWeight, minSurroundConeWeightRelativity);

    [centerLineWeightingFunctions, surroundLineWeightingFunctions] = RGCMosaicConstructor.visualize.conePoolingWeightsMap(...
            pdfExportSubDir, figNo, figPos, theMRGCMosaic, theTargetRGCindex, ...
            spatialSupportCenterDegs, ...
            spatialSupportTickSeparationArcMin, ...
            optimizationResultsStruct.theFinalPooledConeIndicesAndWeights, ...
            minCenterConeWeight, ...
            minSurroundConeWeight, ...
            minSurroundConeWeightRelativity, ...
            visualizedConeWeights, ...
            'contourGenerationMethod', contourGenerationMethod, ...
            'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
            'scaleBarDegs', scaleBarDegs, ...
            'plotTitle', plotTitle, ...
            'axesToRenderIn', axCenterSurroundCompositeConePoolingMap);
 

	% The cell's center and surround pooled cone line weighting functions
	figNo = 11; figPos = [850 460];

	% Visualize the horizontal or vertical line weighting function of the center and surround cone pooling subregion
	whichMeridian = 'horizontal';		% choose between {'horizontal', 'vertical'}
	RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
        pdfExportSubDir, figNo, figPos, ...
		spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
		centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
        'axesToRenderIn', axHorizontalLineWeightingFunctions);

    % The cell's center and surround pooled cone line weighting functions
    figNo = 12; figPos = [1050 460];

    % Visualize the horizontal or vertical line weighting function of the center and surround cone pooling subregion
    whichMeridian = 'vertical';       % choose between {'horizontal', 'vertical'}
    RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
        pdfExportSubDir, figNo, figPos, ...
        spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
        centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
        'axesToRenderIn', axVerticalLineWeightingFunctions);

    if (summaryFigureInsteadOfSeparateFigures)
        switch (targetRFcenterDominantConeType)
            case cMosaic.LCONE_ID
                targetConeDominance = 'LconeDominated';
            case cMosaic.MCONE_ID
                targetConeDominance = 'MconeDominated';
            case cMosaic.SCONE_ID
                targetConeDominance = 'SconeDominated';
            otherwise
                error('Invalid targetRFcenterDominantConeType: %d', targetRFcenterDominantConeType)
        end

        % OLD WAY
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, theSummaryPDFfileName);

        thePDFfileName = fullfile(pdfExportSubDir, theSummaryPDFfileName);
        NicePlot.exportFigToPDF(thePDFfileName, hFigSummary, 300);
    end

end


function renderAchievedOptimizationPerformance(pdfExportSubDir, figNo, figPos, targetVisualSTFparams, ...
        fittedDoGModelParams, fittedDoGModelResidual, totalResidual,...
        fixedH1CellIndex, targetRFcenterConesNum, targetRFcenterDominantConeType, varargin)
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noTitle', false, @islogical);
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noTitle = p.Results.noTitle;
    axesToRenderIn = p.Results.axesToRenderIn;

    targetRsRcRatio = targetVisualSTFparams.surroundToCenterRcRatio;
    targetSCintSensRatio = targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio;
    targetKsKcRatio = targetSCintSensRatio / (targetRsRcRatio^2);

    idx = find(strcmp(fittedDoGModelParams.names, 'RsToRc'));
    achievedRsRcRatio = fittedDoGModelParams.finalValues(idx);

    idx = find(strcmp(fittedDoGModelParams.names, 'intStoCsens'));
    achievedSCintSensRatio = fittedDoGModelParams.finalValues(idx);
    achievedKsKcRatio = achievedSCintSensRatio  / achievedRsRcRatio^2 ;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');

    if (isempty(axesToRenderIn))
        % Initialize figure
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    % Plot
    bar(ax, 1, ((achievedRsRcRatio/targetRsRcRatio)-1), 0.5, 'BaseValue', 0.0);
    hold(ax, 'on');
    bar(ax, 2, ((achievedSCintSensRatio/targetSCintSensRatio)-1), 0.5, 'BaseValue', 0.0);
    bar(ax, 3, fittedDoGModelResidual, 0.5, 'BaseValue', 0.0);
    grid(ax, 'on');
    box(ax, 'off');
    
    if (~noYLabel)
        ylabel(ax,'residuals (achieved/target-1)');
    end
    xlabel(ax, '');

    axis (ax, 'square');
    XLims = [0.75 3.25]; YLims = [-0.2 0.2];
    set(ax, 'XLim', XLims, 'YLim', YLims, ...
            'XTick', 0.5:0.5:3.5, 'XTickLabel', {'', 'R_s/R_c', '', 'intS_s/intS_c', '', 'fit', ''}, ...
            'YTick', -0.5:0.05:0.5);

    % Plot title
    if (~noTitle)
        switch (targetRFcenterDominantConeType)
            case cMosaic.LCONE_ID
                title(ax, sprintf('H1-%s, %d center cones, L-cone dominance\nnorm (residuals) :%2.4f',fixedH1CellIndex, targetRFcenterConesNum, totalResidual));
            case cMosaic.MCONE_ID
                title(ax, sprintf('H1-%s, %d center cones, M-cone dominance\nnorm (residuals) :%2.4f',fixedH1CellIndex, targetRFcenterConesNum, totalResidual));
            case cMosaic.SCONE_ID
                title(ax, sprintf('H1-%s, %d center cones, S-cone dominance\nnorm (residuals) :%2.4f',fixedH1CellIndex, targetRFcenterConesNum, totalResidual));
            otherwise
                error('Incorrect cone type: %d', targetRFcenterDominantConeType)
        end
    end

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, [], []);

    if (isempty(axesToRenderIn))
        % Export figure
        % OLD WAY
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'achievedPerformance.pdf');

        thePDFfileName = fullfile(pdfExportSubDir, 'achievedPerformance.pdf');

        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end




function renderCorrespondenceBetweenFittedH1paramsAndPackerDaceyData(pdfExportSubDir, figNo, figPos,  ...
 		surroundPoolingModelVariables, surroundPoolingModel, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');

    if (isempty(axesToRenderIn))
	   hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    idx = find(strcmp(surroundPoolingModelVariables.names,  'VnVwRatio'));
    fittedModel.NWvolumeRatio = surroundPoolingModelVariables.finalValues(idx);
    
    idx = find(strcmp(surroundPoolingModelVariables.names,  'RnRwRatio'));
    fittedModel.RnarrowToRwideRatio = surroundPoolingModelVariables.finalValues(idx);

    % Outline the tolerance region
    hold(ax, 'on');
    if (contains(surroundPoolingModel.name, 'Fixed'))
        xx = surroundPoolingModel.H1parameterTolerances.RnarrowToRwideRatio*[-1 -1 1 1 -1];
        yy = surroundPoolingModel.H1parameterTolerances.VnarrowToVwideRatio*[-1 1 1 -1 -1]; 
        zz = -10*eps*ones(size(yy));
        for h1Index = 1:numel(RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.RnarrowToRwideRatios)
            xxx = xx + RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.RnarrowToRwideRatios(h1Index);
            yyy = yy + RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.VnarrowToVwideRatios(h1Index);
            patch(ax,xxx,yyy,zz,'FaceColor',[1 1 0.5],'EdgeColor', [1 1 0]*0.5, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
        end
    else
        xStart = surroundPoolingModel.H1parameterTolerances.RnarrowToRwideRatio(1);
        xEnd = surroundPoolingModel.H1parameterTolerances.RnarrowToRwideRatio(2);
        yStart = surroundPoolingModel.H1parameterTolerances.VnarrowToVwideRatio(1);
        yEnd = surroundPoolingModel.H1parameterTolerances.VnarrowToVwideRatio(2);
        xxx = [xStart xStart xEnd xEnd xStart];
        yyy = [yStart yEnd yEnd yStart yStart];
        zz = -10*eps*ones(size(yyy));
        patch(ax,xxx,yyy,zz,'FaceColor',[1 1 0.5],'EdgeColor', [1 1 0]*0.5, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
    end

    % Plot Packer&Dacey H1 data
    plot(ax, RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.RnarrowToRwideRatios,...
             RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.VnarrowToVwideRatios, ...
             'kh', ...
             'MarkerSize', ff.markerSize, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', ff.lineWidth);

    % Plot fitted cell data
    scatter(ax, fittedModel.RnarrowToRwideRatio, fittedModel.NWvolumeRatio, (ff.markerSize)^2,'o', ...
             'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], ...
             'LineWidth', ff.lineWidth*1.5);

	%legend(ax, {surroundPoolingModel.params.PackerDaceyParamsForH1CellIndices.source, 'fitted H1 model'}, ...
    %        'FontSize', ff.legendFontSize, 'box', 'off', 'Color', [0.8 0.8 0.7], ...
    %        'Location', 'NorthOutside');

    xlabel(ax, 'Rnarrow / Rwide');
    ylabel(ax, 'Vnarrow / Vwide');
   

    axis(ax, 'square');
    XLims = [0.0 0.4];
    YLims = [0.0 1.4];
    set(ax, 'XLim', XLims, 'YLim', YLims, ...
    	'YTick', 0:0.2:1.2, ...
    	'XTick', 0:0.05:0.4, 'XTickLabel', {'0', '', '.10', '', '.20', '', '.30', '', '.40'});


   	% Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    if (isempty(axesToRenderIn))
        % Export figure
        % OLD WAY
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'fitted_H1model.pdf');

        thePDFfileName = fullfile(pdfExportSubDir, 'fitted_H1model.pdf');
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end

function renderModelParamsAndRanges(pdfExportSubDir, figNo, figPos, ...
		modelParams, modelName, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    ff.axisFontSize = 18;
    ff.axisTickAngle = 90;

    if (isempty(axesToRenderIn))
    	hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    % Plot
    RGCMosaicConstructor.visualize.fittedModelParams(ax, modelParams, modelName);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    if (isempty(axesToRenderIn))
        % Export figure
        % OLD WAY
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('fitted_%s.pdf', modelName));

        thePDFfileName = fullfile(pdfExportSubDir, sprintf('fitted_%s.pdf', modelName));
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end

end

function renderSTF(pdfExportSubDir, figNo, figPos, sf, STF, sfHiRes, compositeSTFHiRes, centerSTFHiRes, surroundSTFHiRes, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure no left axis label');

    if (isempty(axesToRenderIn))
    	hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    maxSTF = max([max(centerSTFHiRes(:)) max(surroundSTFHiRes(:)) max(compositeSTFHiRes(:)) max(STF(:))]);
    % Plot
    plot(ax, sfHiRes, centerSTFHiRes/maxSTF, 'k-', 'Color', [0.9 0.2 0.5], 'LineWidth', ff.lineWidth*2);
    hold(ax, 'on');
    plot(ax, sfHiRes, surroundSTFHiRes/maxSTF, 'k--', 'Color', [0.3 0.5 0.65], 'LineWidth', ff.lineWidth*2);
    plot(ax, sfHiRes, compositeSTFHiRes/maxSTF, 'k-', 'Color', [0.75 0.75 0.75]*0.75, 'LineWidth', ff.lineWidth*2);

    scatter(ax, sf, STF/maxSTF, (ff.markerSize+2)^2,'o', ...
             'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', [0.75 0.75 0.75]*0.5, ...
             'LineWidth', ff.lineWidth);
    axis(ax, 'square');
    XLims = [0.01 100];
    YLims = [0 1.0];
    set(ax, 'XScale', 'log', 'XLim', XLims, ...
        'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], ...
        'XTickLabel', {'.01', '.03', '.1', '.3', '1', '3', '10', '30', '100'});
    set(ax, 'YLim', YLims, 'YTick', 0:0.1:1.0, 'YTickLabel', {});
    grid(ax, 'on'); box(ax, 'off');
    xlabel(ax, 'spatial frequency (c/deg)');
    ylabel(ax, 'response modulation');
    xtickangle(ax, 0)

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    if (isempty(axesToRenderIn))
        % Export figure
        % OLD WAY
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'fitted_optimizedSTF.pdf');

        thePDFfileName = fullfile(pdfExportSubDir, 'fitted_optimizedSTF.pdf');
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end