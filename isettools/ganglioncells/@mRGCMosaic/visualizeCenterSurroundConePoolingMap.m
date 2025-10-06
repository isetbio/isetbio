function [hFig, ax, centerLineWeightingFunctions, surroundLineWeightingFunctions] = visualizeCenterSurroundConePoolingMap(obj, theRGCindex, varargin)
    
    % Parse input
    p = inputParser;
    p.addParameter('pooledConeIndicesAndWeights', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('visualizedConeWeights', 'surround-alone', @(x)(ismember(x, {'center-surround', 'surround-alone'})));
    p.addParameter('minConeWeightForVisualizingRFcenterPooling', 0.01, @isscalar);
    p.addParameter('minConeWeightForVisualizingRFsurroundPooling', 0.001, @isscalar);
    p.addParameter('minSurroundConeWeightRelativity', 'center', @(x)(ismember(x, {'center', 'surround'})));
    p.addParameter('spatialSupportCenterDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('spatialSupportTickSeparationArcMin', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('withLineWeightingFunctions', false, @islogical);
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('scaleBarDegs', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('doNotLabelScaleBar', false, @islogical);
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConeApertureImage', @(x)(ismember(x, mRGCMosaic.validRFsubregionContourGenerationMethods)));
    p.addParameter('maxNumberOfConesOutsideContour', 1, @isscalar);
    p.addParameter('gridless', false, @islogical);
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.addParameter('figNo', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('figPos', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('withCustomFigureFormat', '', @(x)(isempty(x)||(ischar(x))));
    p.addParameter('pdfExportSubDir', '', @ischar);
    p.addParameter('exportToFigurePDFsDirWithPDFFileName', '', @(x)(isempty(x)||(ischar(x))));
    p.parse(varargin{:});

    pooledConeIndicesAndWeights = p.Results.pooledConeIndicesAndWeights;
    visualizedConeWeights = p.Results.visualizedConeWeights;
    minConeWeightForVisualizingRFcenterPooling = p.Results.minConeWeightForVisualizingRFcenterPooling;
    minConeWeightForVisualizingRFsurroundPooling = p.Results.minConeWeightForVisualizingRFsurroundPooling;
    minSurroundConeWeightRelativity = p.Results.minSurroundConeWeightRelativity;

    spatialSupportCenterDegs = p.Results.spatialSupportCenterDegs;
    spatialSupportTickSeparationArcMin = p.Results.spatialSupportTickSeparationArcMin;

    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    scaleBarDegs = p.Results.scaleBarDegs;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;

    plotLineWeightingFunctions = p.Results.withLineWeightingFunctions;
    gridless = p.Results.gridless;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    plotTitle = p.Results.plotTitle;
    axesToRenderIn = p.Results.axesToRenderIn;
    customFigureFormat = p.Results.withCustomFigureFormat;
    contourGenerationMethod = p.Results.contourGenerationMethod;
    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;

    figNo = p.Results.figNo;
    figPos = p.Results.figPos;
    exportToFigurePDFsDirWithPDFFileName = p.Results.exportToFigurePDFsDirWithPDFFileName;
    pdfExportSubDir = p.Results.pdfExportSubDir;



    [scaleBarDegsDefault, ...
     scaleBarMicronsDefault, ...
	 spatialSupportTickSeparationArcMinDefault, ...
	 spatialSupportCenterDegsDefault, ...
     domainVisualizationLimitsDefault, ...
     domainVisualizationTicksDefault, ...
     domainVisualizationLimitsSingleRFDefault, ...
     domainVisualizationTicksSingleRFDefault] = ...
		 	RGCMosaicAnalyzer.visualize.generateLimits(obj, obj.rgcRFpositionsDegs(theRGCindex,:), ...
            'spatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin);


    if (isempty(spatialSupportCenterDegs))
        spatialSupportCenterDegs = spatialSupportCenterDegsDefault;
    end

    if (isempty(domainVisualizationLimits))
        domainVisualizationLimits = domainVisualizationLimitsSingleRFDefault;
    end

    if (isempty(domainVisualizationTicks))
        domainVisualizationTicks = domainVisualizationTicksSingleRFDefault;
    end

   
    if (isempty(spatialSupportTickSeparationArcMin))
        spatialSupportTickSeparationArcMin = spatialSupportTickSeparationArcMinDefault;
    else
        domainVisualizationTicks
        domainVisualizationLimits
    end


    

    if (isempty(scaleBarDegs))
        scaleBarDegs = scaleBarDegsDefault;
    end



    if (plotLineWeightingFunctions)
        ff = PublicationReadyPlotLib.figureComponents('2x2 standard figure');
    else
        if (~isempty(customFigureFormat))
            ff = PublicationReadyPlotLib.figureComponents(customFigureFormat);
        else
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
        end
    end


    if (isempty(axesToRenderIn))
        if (isempty(figNo))
            figNo = 1000;
        end
	    hFig = figure(figNo); clf;
        figNo = figNo + 1;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
        hFig = [];
    end


    if (plotLineWeightingFunctions)
        axHorizontalProfile = theAxes{2,1};
        axVerticalProfile = theAxes{1,2};
        delete(theAxes{2,2});
    end

    if (isempty(spatialSupportCenterDegs))
        spatialSupportCenterDegs = obj.rgcRFpositionsDegs(theRGCindex,:);
    end


    if (isempty(pooledConeIndicesAndWeights))
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        pooledConeIndicesAndWeights.centerConeIndices = find(connectivityVector > 0);
        pooledConeIndicesAndWeights.centerConeWeights = connectivityVector(pooledConeIndicesAndWeights.centerConeIndices);
    
        connectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
        pooledConeIndicesAndWeights.surroundConeIndices = find(connectivityVector > 0);
        pooledConeIndicesAndWeights.surroundConeWeights = connectivityVector(pooledConeIndicesAndWeights.surroundConeIndices);
    end



    % Surround line weighting functions (using all cones, not just those above minConeWeightForVisualizingRFsurroundPooling)
    surroundLineWeightingFunctions = mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(ax, ...
        obj.inputConeMosaic, [], ...
        pooledConeIndicesAndWeights.surroundConeIndices, ...
        pooledConeIndicesAndWeights.surroundConeWeights, ...
        1.0, spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
        scaleBarDegs, gridless, noXLabel, noYLabel, ...
        'computeLineWeightingFunctionsAndReturn', true, ...
        'doNotLabelScaleBar', doNotLabelScaleBar);

    % Center line weighting functions  (using all cones, not just those above minConeWeightForVisualizingRFcenterPooling)
    centerLineWeightingFunctions = mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(ax, ...
        obj.inputConeMosaic, [], ...
        pooledConeIndicesAndWeights.centerConeIndices, ...
        pooledConeIndicesAndWeights.centerConeWeights, ...
        1.0, spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
        scaleBarDegs, gridless, noXLabel, noYLabel, ...
        'computeLineWeightingFunctionsAndReturn', true, ...
        'doNotLabelScaleBar', doNotLabelScaleBar);


    % Generate visualization cache for this single cell 
    centerSubregionContourSamples = 32;
    contourGenerationMethod = 'ellipseFitToPooledConePositions';
    minConeWeightIncluded = minConeWeightForVisualizingRFcenterPooling;
    obj.generateVisualizationCache([],[], centerSubregionContourSamples, ...
        contourGenerationMethod, theRGCindex, minConeWeightIncluded, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour);


    % XY lims
    spatialSupportRangeArcMin = 60*(max(centerLineWeightingFunctions.xProfile.spatialSupportDegs(:)) - min(centerLineWeightingFunctions.xProfile.spatialSupportDegs(:)));

    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;
    spatialSupportXYDegs(:,1) = spatialSupportCenterDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = spatialSupportCenterDegs(2) + spatialSupportDegs;
    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = spatialSupportCenterDegs(2) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];


    if (isempty(domainVisualizationLimits))
        domainVisualizationLimits = [XLims(1) XLims(2) YLims(1) YLims(2)];
    end

    if (isempty(domainVisualizationTicks))
        domainVisualizationTicks = struct(...
            'x', spatialSupportCenterDegs(1) + (-3:3)*spatialSupportTickSeparationArcMin/60, ...
            'y', spatialSupportCenterDegs(2) + (-3:3)*spatialSupportTickSeparationArcMin/60);
    end


    switch (visualizedConeWeights)
        case 'center-surround'
            % Find cones whose center-surround weight > 0 (center cones in electrophysiological assessment)
            differentialCenterWeights = pooledConeIndicesAndWeights.centerConeWeights;
    
            parfor iCenterConeIndex = 1:numel(pooledConeIndicesAndWeights.centerConeIndices)
                % Find the index of this center cone in the pool of centerConeIndices
                iCenterConeIndexInSurroundPool = ...
                    find(pooledConeIndicesAndWeights.surroundConeIndices == pooledConeIndicesAndWeights.centerConeIndices(iCenterConeIndex));
                if (~isempty(iCenterConeIndexInSurroundPool))
                    differentialCenterWeights(iCenterConeIndex) = pooledConeIndicesAndWeights.centerConeWeights(iCenterConeIndex) - ...
                        pooledConeIndicesAndWeights.surroundConeWeights(iCenterConeIndexInSurroundPool);
                end
            end
    
            % Only identify center cones whose differential weight > minConeWeightForVisualizingRFcenterPooling
            idx = find(differentialCenterWeights > max(differentialCenterWeights) * minConeWeightForVisualizingRFcenterPooling);
            centerAppearingConeIndicesInElectroPhys = pooledConeIndicesAndWeightscenterConeIndices(idx);
   
    
            % Now find cones whose surround-center weight > (surround cones in electrophysiological assessment)
            differentialSurroundWeights = pooledConeIndicesAndWeights.surroundConeWeights;
    
            parfor iSurroundConeIndex = 1:numel(pooledConeIndicesAndWeights.surroundConeIndices)
                % Find the index of this surround cone in the pool of surroundConeIndices
                iSurroundConeIndexInCenterPool = ...
                    find(pooledConeIndicesAndWeights.centerConeIndices == pooledConeIndicesAndWeights.surroundConeIndices(iSurroundConeIndex));
                if (~isempty(iSurroundConeIndexInCenterPool))
                    differentialSurroundWeights(iSurroundConeIndex) = pooledConeIndicesAndWeights.surroundConeWeights(iSurroundConeIndex) - ...
                        pooledConeIndicesAndWeights.centerConeWeights(iSurroundConeIndexInCenterPool);
                end
            end
    
            switch (minSurroundConeWeightRelativity)
                case 'center'
                     maxValue = max(differentialCenterWeights);
                case 'surround'
                     maxValue = max(differentialSurroundWeights);
                otherwise
                     error('Unknown minSurroundConeWeightRelativity: ''%s''.', minSurroundConeWeightRelativity);
            end
    
            % Only identify surround weights whose weight minConeWeightForVisualizingRFsurroundPooling
            idx = find(differentialSurroundWeights > maxValue * minConeWeightForVisualizingRFsurroundPooling);
            surroundAppearingConeIndicesInElectroPhys = pooledConeIndicesAndWeights.surroundConeIndices(idx);
    
            % Visualize both sets of identified cone indices
            identifiedInputConeIndices = cat(1, centerAppearingConeIndicesInElectroPhys(:), surroundAppearingConeIndicesInElectroPhys(:));

    case 'surround-alone'
        % Component weights: surround alone
        % Identify cones pooled by the RF surround
        switch (minSurroundConeWeightRelativity)
            case 'center'
                 maxValue = max(pooledConeIndicesAndWeights.centerConeWeights);
            case 'surround'
                 maxValue = max(pooledConeIndicesAndWeights.surroundConeWeights);
            otherwise
                    error('Unknown minSurroundConeWeightRelativity: ''%s''.', minSurroundConeWeightRelativity)
        end

        % The identified cones
        if (isinf(minConeWeightForVisualizingRFsurroundPooling))
            % Identify cones pooled by the RF center
            idx = find(pooledConeIndicesAndWeights.centerConeWeights >=  maxValue * minConeWeightForVisualizingRFcenterPooling);
            fprintf('Infinite minConeWeightForVisualizingRFsurroundPooling. Only visualizing surround cones\n');
            identifiedInputConeIndices = pooledConeIndicesAndWeights.centerConeIndices(idx);
        else
            idx = find(pooledConeIndicesAndWeights.surroundConeWeights >= maxValue * minConeWeightForVisualizingRFsurroundPooling);
            fprintf('Will visualize %d surround cones (at the %2.3f x ''%s'' level)\n', numel(idx), minConeWeightForVisualizingRFsurroundPooling, minSurroundConeWeightRelativity);
            identifiedInputConeIndices = pooledConeIndicesAndWeights.surroundConeIndices(idx); 
        end

    otherwise
        error('Unknown visualizedConeWeights: ''%s''.', visualizedConeWeights);
    end % switch


    % Plot cones with nice colors
    obj.inputConeMosaic.lConeColor = RGCMosaicConstructor.constants.LcenterColor;
    obj.inputConeMosaic.mConeColor = RGCMosaicConstructor.constants.McenterColor;
    obj.inputConeMosaic.sConeColor = RGCMosaicConstructor.constants.ScenterColor;


    obj.visualize(...
        'axesHandle', ax, ...
        'scaleBarDegs', scaleBarDegs, ...
        'doNotLabelScaleBar', doNotLabelScaleBar, ...
        'visualizedRGCindices', theRGCindex, ...
        'pooledConesLineWidth', [], ...
        'plottedRFoutlineLineWidth', 4.0, ...
        'plottedRFoutlineFaceAlpha', 0.5, ...
        'identifyPooledCones', true, ...
        'identifyInputCones', true, ...
        'inputConesAlpha', 1.0, ...
        'identifiedConeApertureThetaSamples', 40, ...
        'identifiedConeAperture', 'lightCollectingArea4sigma', ...
        'identifiedInputConeIndices', identifiedInputConeIndices, ...                 % which cones to identity with their type
        'identifiedInputConeIndicesContour', true, ...
        'minConeWeightVisualized', minConeWeightForVisualizingRFcenterPooling, ...    % where to draw the center profile
        'contourGenerationMethod', contourGenerationMethod, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
        'centerSubregionContourSamples', centerSubregionContourSamples, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'noXLabel', plotLineWeightingFunctions, ...
        'plotTitle', plotTitle ...
        );


    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    if (plotLineWeightingFunctions)
        % The horizontal profile
        whichMeridian = 'horizontal';		% choose between {'horizontal', 'vertical'}
	    RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
            pdfExportSubDir, figNo, figPos, ...
		    spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
		    centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
            'domainVisualizationLimits', domainVisualizationLimits(1:2), ...
            'domainVisualizationTicks', domainVisualizationTicks.x, ...
            'axesToRenderIn', axHorizontalProfile);
        PublicationReadyPlotLib.applyFormat(axHorizontalProfile,ff);
        figNo = figNo + 1;

        % The vertical profile
        whichMeridian = 'vertical';		% choose between {'horizontal', 'vertical'}
	    RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
            pdfExportSubDir, figNo, figPos, ...
		    spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
		    centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
            'domainVisualizationLimits', domainVisualizationLimits(3:4), ...
            'domainVisualizationTicks', domainVisualizationTicks.y, ...
            'axesToRenderIn', axVerticalProfile, ...
            'flipXYaxes', true);
        PublicationReadyPlotLib.applyFormat(axVerticalProfile,ff);
        figNo = figNo + 1;
    end % if (plotLineWeightingFunctions)

    
    if (isempty(axesToRenderIn))
        % Export figure

        if (~isempty(exportToFigurePDFsDirWithPDFFileName))
            p = getpref('isetbio');
            pdfExportRootDir = fullfile(p.rgcResources.figurePDFsDir);
            theVisualizationPDFfilename = fullfile(pdfExportSubDir, exportToFigurePDFsDirWithPDFFileName);
    
            % Generate the path if we need to
            RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);
    
            thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
            NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);
        else
            thePDFFileName = sprintf('conePoolingWeightsMap_centerConeThreshold_%2.2f_surroundConeThreshold_%2.3f.pdf', ...
                minConeWeightForVisualizingRFcenterPooling, ...
                minConeWeightForVisualizingRFsurroundPooling);
    
            % OLD WAY
            %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFFileName);
    
            thePDFfileName = fullfile(pdfExportSubDir, thePDFFileName);
            NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
        end

    end

end
