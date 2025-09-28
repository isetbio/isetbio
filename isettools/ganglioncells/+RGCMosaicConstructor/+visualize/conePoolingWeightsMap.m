function [centerLineWeightingFunctions, surroundLineWeightingFunctions] = ...
	conePoolingWeightsMap(pdfExportSubDir, figNo, figPos, theMRGCMosaic, theTargetRGCindex, ...
		spatialSupportCenterDegs, spatialSupportTickSeparationArcMin,  ...
		theOptimizedPooledConeIndicesAndWeights, ...
		minConeWeightForVisualizingRFcenterPooling, ...
        minConeWeightForVisualizingRFsurroundPooling, ...
        minSurroundConeWeightRelativity, ...
        visualizedConeWeights, ...
        varargin)
	
    
    % Parse input
    p = inputParser;
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
    p.addParameter('exportToFigurePDFsDirWithPDFFileName', '', @(x)(isempty(x)||(ischar(x))));
    p.parse(varargin{:});

    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    scaleBarDegs = p.Results.scaleBarDegs;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;

    gridless = p.Results.gridless;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    plotTitle = p.Results.plotTitle;
    axesToRenderIn = p.Results.axesToRenderIn;
    contourGenerationMethod = p.Results.contourGenerationMethod;
    maxNumberOfConesOutsideContour = p.Results.maxNumberOfConesOutsideContour;
    exportToFigurePDFsDirWithPDFFileName = p.Results.exportToFigurePDFsDirWithPDFFileName;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');

    if (isempty(axesToRenderIn))
	    hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end


    % Surround line weighting functions (using all cones, not just those above minConeWeightForVisualizingRFsurroundPooling)
    surroundLineWeightingFunctions = mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(ax, ...
        theMRGCMosaic.inputConeMosaic, [], ...
        theOptimizedPooledConeIndicesAndWeights.surroundConeIndices, ...
        theOptimizedPooledConeIndicesAndWeights.surroundConeWeights, ...
        1.0, spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
        scaleBarDegs, gridless, noXLabel, noYLabel, ...
        'computeLineWeightingFunctionsAndReturn', true, ...
        'doNotLabelScaleBar', doNotLabelScaleBar);

    % Center line weighting functions  (using all cones, not just those above minConeWeightForVisualizingRFcenterPooling)
    centerLineWeightingFunctions = mRGCMosaic.renderInputConeMosaicSubregionPoolingMap(ax, ...
        theMRGCMosaic.inputConeMosaic, [], ...
        theOptimizedPooledConeIndicesAndWeights.centerConeIndices, ...
        theOptimizedPooledConeIndicesAndWeights.centerConeWeights, ...
        1.0, spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
        scaleBarDegs, gridless, noXLabel, noYLabel, ...
        'computeLineWeightingFunctionsAndReturn', true, ...
        'doNotLabelScaleBar', doNotLabelScaleBar);

    % Generate visualization cache for this single cell 
    centerSubregionContourSamples = 32;
    contourGenerationMethod = 'ellipseFitToPooledConePositions';
    minConeWeightIncluded = minConeWeightForVisualizingRFcenterPooling;
    theMRGCMosaic.generateVisualizationCache([],[], centerSubregionContourSamples, ...
        contourGenerationMethod, theTargetRGCindex, minConeWeightIncluded, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour);

    % The cell's RF center contour
    theContourData = theMRGCMosaic.visualizationCache.rfCenterContourData{1};
    % XY lims
    spatialSupportRangeArcMin = 3 * spatialSupportTickSeparationArcMin;
    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;
    spatialSupportXYDegs(:,1) = spatialSupportCenterDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = spatialSupportCenterDegs(2) + spatialSupportDegs;
    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = spatialSupportCenterDegs(2) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];

    if (isempty(domainVisualizationLimits))
        domainVisualizationLimits = [XLims(1) XLims(2) YLims(1) YLims(2) ];
    end

    if (isempty(domainVisualizationTicks))
        domainVisualizationTicks = struct(...
            'x', spatialSupportCenterDegs(1) + (-3:3)*spatialSupportTickSeparationArcMin/60, ...
            'y', spatialSupportCenterDegs(2) + (-3:3)*spatialSupportTickSeparationArcMin/60);
    end


    switch (visualizedConeWeights)
        case 'center-surround'
            % Find cones whose center-surround weight > 0 (center cones in electrophysiological assessment)
            differentialCenterWeights = theOptimizedPooledConeIndicesAndWeights.centerConeWeights;
    
            parfor iCenterConeIndex = 1:numel(theOptimizedPooledConeIndicesAndWeights.centerConeIndices)
                % Find the index of this center cone in the pool of theOptimizedPooledConeIndicesAndWeights.surroundConeIndices
                iCenterConeIndexInSurroundPool = ...
                    find(theOptimizedPooledConeIndicesAndWeights.surroundConeIndices == theOptimizedPooledConeIndicesAndWeights.centerConeIndices(iCenterConeIndex));
                if (~isempty(iCenterConeIndexInSurroundPool))
                    differentialCenterWeights(iCenterConeIndex) = theOptimizedPooledConeIndicesAndWeights.centerConeWeights(iCenterConeIndex) - ...
                        theOptimizedPooledConeIndicesAndWeights.surroundConeWeights(iCenterConeIndexInSurroundPool);
                end
            end
    
            % Only identify center cones whose differential weight > minConeWeightForVisualizingRFcenterPooling
            idx = find(differentialCenterWeights > max(differentialCenterWeights) * minConeWeightForVisualizingRFcenterPooling);
            centerAppearingConeIndicesInElectroPhys = theOptimizedPooledConeIndicesAndWeights.centerConeIndices(idx);
    
            
    
            % Now find cones whose surround-center weight > (surround cones in electrophysiological assessment)
            differentialSurroundWeights = theOptimizedPooledConeIndicesAndWeights.surroundConeWeights;
    
            parfor iSurroundConeIndex = 1:numel(theOptimizedPooledConeIndicesAndWeights.surroundConeIndices)
                % Find the index of this surround cone in the pool of theOptimizedPooledConeIndicesAndWeights.centerConeIndices
                iSurroundConeIndexInCenterPool = ...
                    find(theOptimizedPooledConeIndicesAndWeights.centerConeIndices == theOptimizedPooledConeIndicesAndWeights.surroundConeIndices(iSurroundConeIndex));
                if (~isempty(iSurroundConeIndexInCenterPool))
                    differentialSurroundWeights(iSurroundConeIndex) = theOptimizedPooledConeIndicesAndWeights.surroundConeWeights(iSurroundConeIndex) - ...
                        theOptimizedPooledConeIndicesAndWeights.centerConeWeights(iSurroundConeIndexInCenterPool);
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
            surroundAppearingConeIndicesInElectroPhys = theOptimizedPooledConeIndicesAndWeights.surroundConeIndices(idx);
    
            % Visualize both sets of identified cone indices
            identifiedInputConeIndices = cat(1, centerAppearingConeIndicesInElectroPhys(:), surroundAppearingConeIndicesInElectroPhys(:));

    case 'surround-alone'
        % Component weights: surround alone
        % Identify cones pooled by the RF surround
        switch (minSurroundConeWeightRelativity)
            case 'center'
                 maxValue = max(theOptimizedPooledConeIndicesAndWeights.centerConeWeights);
            case 'surround'
                 maxValue = max(theOptimizedPooledConeIndicesAndWeights.surroundConeWeights);
            otherwise
                    error('Unknown minSurroundConeWeightRelativity: ''%s''.', minSurroundConeWeightRelativity)
        end

        % The identified cones
        if (isinf(minConeWeightForVisualizingRFsurroundPooling))
            % Identify cones pooled by the RF center
            idx = find(theOptimizedPooledConeIndicesAndWeights.centerConeWeights >=  maxValue * minConeWeightForVisualizingRFcenterPooling);
            fprintf('Infinite minConeWeightForVisualizingRFsurroundPooling. Only visualizing surround cones\n');
            identifiedInputConeIndices = theOptimizedPooledConeIndicesAndWeights.centerConeIndices(idx);
        else
            idx = find(theOptimizedPooledConeIndicesAndWeights.surroundConeWeights >= maxValue * minConeWeightForVisualizingRFsurroundPooling);
            fprintf('Will visualize %d surround cones (at the %2.3f x ''%s'' level)\n', numel(idx), minConeWeightForVisualizingRFsurroundPooling, minSurroundConeWeightRelativity);
            identifiedInputConeIndices = theOptimizedPooledConeIndicesAndWeights.surroundConeIndices(idx); 
        end

    otherwise
        error('Unknown visualizedConeWeights: ''%s''.', visualizedConeWeights);
    end % switch

    % Plot cones with nice colors, consistent with the STFs
    theMRGCMosaic.inputConeMosaic.lConeColor = RGCMosaicConstructor.constants.LcenterColor;
    theMRGCMosaic.inputConeMosaic.mConeColor = RGCMosaicConstructor.constants.McenterColor;
    theMRGCMosaic.inputConeMosaic.sConeColor = RGCMosaicConstructor.constants.ScenterColor;

    theMRGCMosaic.visualize(...
        'axesHandle', ax, ...
        'scaleBarDegs', scaleBarDegs, ...
        'doNotLabelScaleBar', doNotLabelScaleBar, ...
        'visualizedRGCindices', theTargetRGCindex, ...
        'pooledConesLineWidth', [], ...
        'plottedRFoutlineLineWidth', 4.0, ...
        'plottedRFoutlineFaceAlpha', 0.5, ...
        'identifyPooledCones', true, ...
        'identifyInputCones', true, ...
        'inputConesAlpha', 1.0, ...
        'identifiedConeApertureThetaSamples', 40, ...
        'identifiedConeAperture', 'lightCollectingArea4sigma', ...
        'identifiedInputConeIndices', identifiedInputConeIndices, ...                 % which cones to identity with their type
        'minConeWeightVisualized', minConeWeightForVisualizingRFcenterPooling, ...    % where to draw the center profile
        'contourGenerationMethod', contourGenerationMethod, ...
        'maxNumberOfConesOutsideContour', maxNumberOfConesOutsideContour, ...
        'centerSubregionContourSamples', centerSubregionContourSamples, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'plotTitle', plotTitle ...
        );
   
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);


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