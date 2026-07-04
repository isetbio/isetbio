function generateRFcenterOverlap(obj, rfCenterOverlapParams, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('visualizeGenerationOfOverlappingRFcenterWeights', false, @islogical);
    p.addParameter('minConeWeightVisualized', mRGCMosaic.minRFcenterConeWeightIncludedToMatchFigure4OfFieldEtAl2010, @isscalar);
    p.addParameter('visualizeMinSensitivityForInclusion', true, @islogical);
    p.addParameter('visualizeSensitivityAtPointOfOverlap', true, @islogical);
    p.addParameter('visualizeGauthierSensitivityAtPointOfOverlap', true, @islogical);

    p.parse(varargin{:});
    visualizeGenerationOfOverlappingRFcenterWeights = p.Results.visualizeGenerationOfOverlappingRFcenterWeights;
    minConeWeightVisualized = p.Results.minConeWeightVisualized;
    visualizeMinSensitivityForInclusion = p.Results.visualizeMinSensitivityForInclusion;
    visualizeSensitivityAtPointOfOverlap = p.Results.visualizeSensitivityAtPointOfOverlap;
    visualizeGauthierSensitivityAtPointOfOverlap = p.Results.visualizeGauthierSensitivityAtPointOfOverlap;

    % Validate rfCenterOverlapParams
    assert(isstruct(rfCenterOverlapParams), 'The rfCenterOverlapParams must be a struct');
    assert(isfield(rfCenterOverlapParams, 'minSensitivityForInclusion'),  'The rfCenterOverlapParams must have a ''minSensitivityForInclusion'' scalar field.');
    assert(isfield(rfCenterOverlapParams, 'maxNumberOfConesOutsideContour'),  'The rfCenterOverlapParams must have a ''maxNumberOfConesOutsideContour'' scalar field.');
    assert(isfield(rfCenterOverlapParams, 'sensitivityAtPointOfOverlap'), 'The rfCenterOverlapParams must have a ''sensitivityAtPointOfOverlap'' scalar field.');
    assert(isfield(rfCenterOverlapParams, 'overlappingWeightsDivergencePattern'), 'The rfCenterOverlapParams must have a ''overlappingWeightsDivergencePattern'' char field.');
    assert(isfield(rfCenterOverlapParams, 'coneTypesIncluded'), 'The rfCenterOverlapParams must have a ''coneTypesIncluded'' field.');
    assert(~isempty(rfCenterOverlapParams.coneTypesIncluded), 'The ''conectableConeTypes'' cannot be empty');
    assert(all(ismember(rfCenterOverlapParams.coneTypesIncluded, [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID])), 'The ''coneTypesIncluded'' are not valid.');
    assert(ismember(rfCenterOverlapParams.overlappingWeightsDivergencePattern, {'isotropic', 'inLineWithExclusiveConnections', 'orthogonalToExclusiveConnections'}), 'The ''overlappingWeightsDivergencePattern'' is not valid.');
    if (obj.rgcRFcentersOverlap)
        fprintf('\nRF centers are already overlapping in this mosaic. Not recomputing center cone weights.\n')
    	return;
    end

	% Save the no-overlap connectivity matrix
	if (isempty(obj.rgcRFcenterConeConnectivityMatrixNoOverlap))
		obj.rgcRFcenterConeConnectivityMatrixNoOverlap = obj.rgcRFcenterConeConnectivityMatrix;
		fprintf('Saved the no-overlap center cone connectivity matrix.\n')
	end


    % Initialize the intermediate cell arrays
    rgcRFcenterOverlappingConeConnectivityIndicesMatrix = cell(1, obj.rgcsNum);
    rgcRFcenterOverlappingConeConnectivityWeightsMatrix = cell(1, obj.rgcsNum);
    exclusivelyConnectedInputConeIndicesNum = zeros(1, obj.rgcsNum);

    % Compute overlapping center cone weights
    fprintf('Computing weights for overlaping centers in %d mRGCs. Please wait ....\n', obj.rgcsNum);

    % Compute super Gaussian exponent for all mRGCs in the mosaic
    gaussianPower = computeSuperGaussianExponent(obj.rgcRFspacingsDegs, obj.rgcRFpositionsDegs, ...
        obj.superGaussianMaxExponent, obj.superGaussianExponentSigmoidSteepnessExponent, ...
        obj.superGaussianExponentSigmoidEccDegs50);
    
    for theRGCindex = 1:obj.rgcsNum

        if (mod(theRGCindex, 100) == 0)
            fprintf('Generating RF overlap for mRGC %d of %d\n', theRGCindex, obj.rgcsNum);
        end

        % Retrieve the indices and weights of the exclusively pooled cones
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrixNoOverlap(:, theRGCindex)));

        indicesOfExclusivelyPooledCones = find(...
            (connectivityVector > 0.001) & ...
            (ismember(obj.inputConeMosaic.coneTypes,rfCenterOverlapParams.coneTypesIncluded)));

        weightsOfExclusivelyPooledCones = connectivityVector(indicesOfExclusivelyPooledCones);
        exclusivelyConnectedInputConeIndicesNum(theRGCindex) = numel(indicesOfExclusivelyPooledCones);

        if (gaussianPower(theRGCindex) >= 0.99*obj.superGaussianMaxExponent)
            % No overlap in RGCs with peak superGaussianExponent
            rgcRFcenterOverlappingConeConnectivityIndicesMatrix{theRGCindex} = [];
            rgcRFcenterOverlappingConeConnectivityWeightsMatrix{theRGCindex} = [];
            continue;
        end


        % Generate the Gaussian ellipsoid interpolant 
        [overlappingRFcenterWeightsInterpolant, xSupport, ySupport, theEllipseRotation] = generateOverlappingWeightsInterpolant(...
            obj.inputConeMosaic.coneRFpositionsDegs, ...
            obj.inputConeMosaic.coneRFspacingsDegs, ...
            obj.inputConeMosaic.coneTypes, ...
            connectivityVector, ...
            indicesOfExclusivelyPooledCones, ...
            rfCenterOverlapParams.sensitivityAtPointOfOverlap, ...
            rfCenterOverlapParams.maxNumberOfConesOutsideContour, ...
            rfCenterOverlapParams.overlappingWeightsDivergencePattern, ...
            rfCenterOverlapParams.coneTypesIncluded, ...
            gaussianPower(theRGCindex), ...
            theRGCindex);

        % Apply Gaussian ellipsoid interpolant to all cones
        overlappingPoolingWeights = overlappingRFcenterWeightsInterpolant(...
            obj.inputConeMosaic.coneRFpositionsDegs(:,1), ...
            obj.inputConeMosaic.coneRFpositionsDegs(:,2));

        % Include cones with weights down to rfCenterOverlapParams.minSensitivityForInclusion
        % Only including cones of the desired type, rfCenterOverlapParams.coneTypesIncluded
        % BUT THIS SHOULD BE chromatic-spatial tradeoff param controlled
        theOverlappingRFcenterPooledConeIndices = find(...
            (overlappingPoolingWeights >= rfCenterOverlapParams.minSensitivityForInclusion) & ...
            (ismember(obj.inputConeMosaic.coneTypes, rfCenterOverlapParams.coneTypesIncluded))...
            );
        theOverlappingRFcenterPoolingWeights = overlappingPoolingWeights(theOverlappingRFcenterPooledConeIndices);

        coneIndicesWithinPointofOverlap = find(theOverlappingRFcenterPoolingWeights > rfCenterOverlapParams.sensitivityAtPointOfOverlap);
        if (numel(coneIndicesWithinPointofOverlap) > exclusivelyConnectedInputConeIndicesNum(theRGCindex))
            fprintf('Extra cones within sensitivity at overlap for RGC %d (%d expected %d). Attenuating their weights.\n ', ...
                theRGCindex, numel(coneIndicesWithinPointofOverlap), exclusivelyConnectedInputConeIndicesNum(theRGCindex));
            
            nCones = exclusivelyConnectedInputConeIndicesNum(theRGCindex)-numel(coneIndicesWithinPointofOverlap)+1;
            [~,idx] = sort(theOverlappingRFcenterPoolingWeights, 'descend');
            theOverlappingRFcenterPoolingWeights = theOverlappingRFcenterPoolingWeights(idx);
            theOverlappingRFcenterPooledConeIndices = theOverlappingRFcenterPooledConeIndices(idx);

            coneIndicesToAdjust = exclusivelyConnectedInputConeIndicesNum(theRGCindex)-(0:nCones);
            gains = 1.001 * rfCenterOverlapParams.sensitivityAtPointOfOverlap ./ theOverlappingRFcenterPoolingWeights(coneIndicesToAdjust);
            theOverlappingRFcenterPoolingWeights(coneIndicesToAdjust) = ...
                theOverlappingRFcenterPoolingWeights(coneIndicesToAdjust) .* gains;

            coneIndicesWithinPointofOverlap = find(theOverlappingRFcenterPoolingWeights > rfCenterOverlapParams.sensitivityAtPointOfOverlap);
            exclusivelyConnectedInputConeIndicesNum(theRGCindex) = numel(coneIndicesWithinPointofOverlap);
        end


       % Sanity check
        if (numel(coneIndicesWithinPointofOverlap) ~= exclusivelyConnectedInputConeIndicesNum(theRGCindex))
            gaussianPower(theRGCindex)
            [min(theOverlappingRFcenterPoolingWeights) max(theOverlappingRFcenterPoolingWeights) rfCenterOverlapParams.sensitivityAtPointOfOverlap]
            [numel(coneIndicesWithinPointofOverlap)  exclusivelyConnectedInputConeIndicesNum(theRGCindex)]
            error('How can this be?')
        end

        % Assign weights to connectivity matrix
        rgcRFcenterOverlappingConeConnectivityIndicesMatrix{theRGCindex} = theOverlappingRFcenterPooledConeIndices;
        rgcRFcenterOverlappingConeConnectivityWeightsMatrix{theRGCindex} = theOverlappingRFcenterPoolingWeights;

        if (visualizeGenerationOfOverlappingRFcenterWeights) 
            tickSeparationArcMin = 10;
            overlapColor = [0.2 0.4 0.7]*1.2;

            hFigSamplerExclusiveWeights = figure(2002); clf;
            set(hFigSamplerExclusiveWeights, 'Name', 'exclusive weights');
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theSamplerAxes = PublicationReadyPlotLib.generatePanelAxes(hFigSamplerExclusiveWeights,ff);

            obj.renderConePoolingRFmap(theRGCindex, ...
                'figureHandle', hFigSamplerExclusiveWeights, ...
                'axesHandle', theSamplerAxes{1,1}, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'forceRegenerateVisualizationCache', true, ...
                'minConeWeightIncluded', minConeWeightVisualized, ...
                'renderSubregionContour', true, ...
                'clearAxes', false, ...
                'renderConeMap', false);
            set(theSamplerAxes{1,1}, 'CLim', [0 1.5])
            colormap(theSamplerAxes{1,1}, brewermap(256, 'greys'));

            XLim = get(theSamplerAxes{1,1}, 'XLim');
            YLim = get(theSamplerAxes{1,1}, 'YLim');
            nSamples = 256;
            xPts = linspace(XLim(1), XLim(2), nSamples);
            yPts = linspace(YLim(1), YLim(2), nSamples);
            [Xq,Yq] = meshgrid(xPts, yPts);
            theSampler = overlappingRFcenterWeightsInterpolant(Xq(:),Yq(:));
            theSampler = theSampler / max(theSampler(:));

            theSamplerMap = reshape(theSampler, nSamples*[1 1]);
            zLevels = [mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections*1.01];
            contour(theSamplerAxes{1,1}, xPts, yPts, theSamplerMap, zLevels, 'LineWidth', 3.0, 'EdgeColor', overlapColor, 'LineStyle', '-');
            contour(theSamplerAxes{1,1}, xPts, yPts, theSamplerMap, zLevels, 'LineWidth', 1.5, 'EdgeColor', [1 1 1], 'LineStyle', '--');
            hold(theSamplerAxes{1,1}, 'on');

            [~,idx] = max(theSamplerMap(:));
            [row, col] = ind2sub(size(theSamplerMap), idx);
            plot(theSamplerAxes{1,1}, xPts, YLim(1) + (YLim(2)-YLim(1)) * 0.02 + (((YLim(2)-YLim(1)) * 0.95)*theSamplerMap(row,:)), '-', ...
                'Color', overlapColor, 'LineWidth', ff.lineWidth);


            ff.box = 'on';
            PublicationReadyPlotLib.applyFormat(theSamplerAxes{1,1},ff);
            pos = get(hFigSamplerExclusiveWeights, 'Position');
            pos(1:2) = [760 700];
            set(hFigSamplerExclusiveWeights, 'Position',pos);

            

            % Overwrite the connectivity matrix so we can visualize the divergence. 
            obj.rgcRFcenterConeConnectivityMatrix(rgcRFcenterOverlappingConeConnectivityIndicesMatrix{theRGCindex}, theRGCindex) = ...
                rgcRFcenterOverlappingConeConnectivityWeightsMatrix {theRGCindex};


            % Plot the overlapping cone pooling (after introducing overlap)
            hFigConePoolingMapAfter = figure(2003); clf;
            set(hFigConePoolingMapAfter, 'Name', 'with overlap');
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theConePoolingMapAfterAxes = PublicationReadyPlotLib.generatePanelAxes(hFigConePoolingMapAfter,ff);


            obj.renderConePoolingRFmap(theRGCindex, ...
                'figureHandle', hFigConePoolingMapAfter, ...
                'axesHandle', theConePoolingMapAfterAxes{1,1}, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'forceRegenerateVisualizationCache', true, ...
                'minConeWeightIncluded', minConeWeightVisualized, ...
                'renderSubregionContour', false);

            hold(theConePoolingMapAfterAxes{1,1}, 'on');
            contour(theConePoolingMapAfterAxes{1,1}, xPts, yPts, theSamplerMap, zLevels, 'LineWidth', 3.0, 'EdgeColor', overlapColor, 'LineStyle', '-');
            contour(theConePoolingMapAfterAxes{1,1}, xPts, yPts, theSamplerMap, zLevels, 'LineWidth', 1.5, 'EdgeColor', [1 1 1], 'LineStyle', '--');

            ff.box = 'on';
            PublicationReadyPlotLib.applyFormat(theConePoolingMapAfterAxes{1,1},ff);
            pos = get(hFigConePoolingMapAfter, 'Position');
            pos(1:2) = [50 125];
            set(hFigConePoolingMapAfter, 'Position',pos);
            

            hFigWeights = figure(2005); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theWeightsAxes = PublicationReadyPlotLib.generatePanelAxes(hFigWeights ,ff);
            xtickIncrement = 1;
            if (numel(theOverlappingRFcenterPoolingWeights) >= 40)
                xtickIncrement = 10;
            elseif (numel(theOverlappingRFcenterPoolingWeights) >= 20)
                xtickIncrement = 5;
            elseif  (numel(theOverlappingRFcenterPoolingWeights) >= 10)
                xtickIncrement = 2;
            end
            
            theWeightBarWidth = 0.3;
            [~,idx] = sort(theOverlappingRFcenterPoolingWeights, 'descend');
            p1 = RGCMosaicAnalyzer.visualize.shadedHistogram(theWeightsAxes{1,1}, ...
                (1:numel(weightsOfExclusivelyPooledCones))-0.5*theWeightBarWidth, weightsOfExclusivelyPooledCones, ...
                theWeightBarWidth, false, [0.6 0.6 0.6], [0 0 0], 0.5, ff.lineWidth, '-');

            legendPlotHandles = [];
            theLegends = {};
            legendPlotHandles(numel(legendPlotHandles)+1) = p1;
            theLegends{numel(theLegends)+1} = 'exclusive';

            hold(theWeightsAxes{1,1}, 'on')
            p2 = RGCMosaicAnalyzer.visualize.shadedHistogram(theWeightsAxes{1,1}, ...
                (1:numel(idx))+0.5*theWeightBarWidth, theOverlappingRFcenterPoolingWeights(idx), ...
                theWeightBarWidth, false, [0.2 0.4 0.7]*1.2, [0.2 0.4 0.7]*0.5, 0.5, ff.lineWidth, '-');
            legendPlotHandles(numel(legendPlotHandles)+1) = p2;
            theLegends{numel(theLegends)+1} = 'overlapping';

            RGCMosaicAnalyzer.visualize.shadedHistogram(theWeightsAxes{1,1}, ...
                (1:numel(weightsOfExclusivelyPooledCones))-0.5*theWeightBarWidth, weightsOfExclusivelyPooledCones, ...
                theWeightBarWidth, false, [0.6 0.6 0.6], [0 0 0], 0.5, ff.lineWidth, '-');

            if (visualizeMinSensitivityForInclusion)
                p3 = plot(theWeightsAxes{1,1}, ...
                    [0 max([6 numel(theOverlappingRFcenterPoolingWeights)])+1] , ...
                    rfCenterOverlapParams.minSensitivityForInclusion*[1 1], ...
                    'g--', 'LineWidth', 1.0);
                legendPlotHandles(numel(legendPlotHandles)+1) = p3;
                theLegends{numel(theLegends)+1} = 'min sens. for inclusion';
            end

            if (visualizeSensitivityAtPointOfOverlap)
                p4 = plot(theWeightsAxes{1,1}, ...
                    [0 max([6 numel(theOverlappingRFcenterPoolingWeights)])+1] , ...
                    rfCenterOverlapParams.sensitivityAtPointOfOverlap*[1 1], ...
                    'r--',  'LineWidth', 1.0);
                legendPlotHandles(numel(legendPlotHandles)+1) = p4;
                theLegends{numel(theLegends)+1} = 'assumed sens. at point of RF overlap';
            end

            if (visualizeGauthierSensitivityAtPointOfOverlap)
                p5 = plot(theWeightsAxes{1,1}, ...
                    [0 max([6 numel(theOverlappingRFcenterPoolingWeights)])+1] , ...
                    mRGCMosaic.sensitivityAtPointOfOverlapFromGauthierRFmappingExperiment*[1 1], ...
                    'k:', 'LineWidth', 1.0);
                legendPlotHandles(numel(legendPlotHandles)+1) = p5;
                theLegends{numel(theLegends)+1} = 'Gauthier sens. at point of RF overlap';
            end

            legend(theWeightsAxes{1,1}, legendPlotHandles, theLegends);

            grid(theWeightsAxes{1,1}, 'on');
            yAxisHandle = get(theWeightsAxes{1,1}, 'YAxis');
            yAxisHandle.TickLabelInterpreter = 'latex';
            xlabel(theWeightsAxes{1,1}, 'cone index');
            ylabel(theWeightsAxes{1,1}, 'pooling weight');
            set(theWeightsAxes{1,1}, 'XLim', [0 max([6 numel(theOverlappingRFcenterPoolingWeights)])+1], 'YLim', [0 1.0], ...
                'YTick', [exp(-4) exp(-2) exp(-1) exp(-0.5) 1], 'YTickLabel', {'$\mathrm e^{-4}$', '$\mathrm e^{-2}$', '$\mathrm e^{-1}$', '$\mathrm e^{-.5}$', '$e^{0}$'}, ...
                'XTick', 0:xtickIncrement:max([6 numel(theOverlappingRFcenterPoolingWeights)]));

            ff.box = 'off';
            PublicationReadyPlotLib.applyFormat(theWeightsAxes{1,1},ff);
            pos = get(hFigWeights, 'Position');
            pos(1:2) = [1350 430];
            set(hFigWeights, 'Position',pos);

            % Export figure
            theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            thePDFfileName = fullfile(theRawFiguresDir, sprintf('ConePoolingMapAfter.pdf'));
            NicePlot.exportFigToPDF(thePDFfileName, hFigConePoolingMapAfter,  300);

            % Export figure
            theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            thePDFfileName = fullfile(theRawFiguresDir, sprintf('GaussianSamplerExclusiveWeights.pdf'));
            NicePlot.exportFigToPDF(thePDFfileName, hFigSamplerExclusiveWeights,  300);
            
            % Export figure
            theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            thePDFfileName = fullfile(theRawFiguresDir, sprintf('Weights.pdf'));
            NicePlot.exportFigToPDF(thePDFfileName, hFigWeights,  300);

            % Pause
            RGCMosaicConstructor.helper.queryUserFor.unpausingExecution(sprintf('Paused by %s', mfilename));

        end % if (visualizeGenerationOfOverlappingRFcenterWeights) 
    end % parfor theRGCindex

    % Initialize the new (overlapping connectivity matrix)
    obj.rgcRFcenterConeConnectivityMatrix = obj.rgcRFcenterConeConnectivityMatrixNoOverlap;

    % Set the exclusivelyConnectedInputConeIndicesNum
    obj.exclusivelyConnectedInputConeIndicesNum = exclusivelyConnectedInputConeIndicesNum;

    % Assign the computed overlapping cone pooling weights
    for theRGCindex = 1:obj.rgcsNum
        theInputConeIndices = rgcRFcenterOverlappingConeConnectivityIndicesMatrix{theRGCindex};
        if (isempty(theInputConeIndices))
            % Must be a cone input numerosity for which we did not allow overlap
            fprintf('No overlapping cone weights for RGC with %d exclusive cone inputs\n', exclusivelyConnectedInputConeIndicesNum(theRGCindex));
        else
            obj.rgcRFcenterConeConnectivityMatrix(:,theRGCindex) = 0*obj.rgcRFcenterConeConnectivityMatrix(:,theRGCindex);
            obj.rgcRFcenterConeConnectivityMatrix(theInputConeIndices, theRGCindex) = rgcRFcenterOverlappingConeConnectivityWeightsMatrix{theRGCindex};
        end
    end

    % Save the rfCenterOverlapParams
    obj.rfCenterOverlapParams = rfCenterOverlapParams;

    % Set the rgcRFcentersOverlap flag to true
    obj.rgcRFcentersOverlap = true;
end


function gaussianPower = computeSuperGaussianExponent(rgcRFspacingsDegs, rgcRFpositionsDegs, ...
    superGaussianMaxExponent, superGaussianExponentSigmoidalVariationExponent, superGaussianExponentSigmoidalVariatioEccDegs50 )
    if (superGaussianMaxExponent < 2)
        error('maxGaussianPower in super Gaussian must be > 2');
    end

    % Convert spacingsDegs to eccentricities empirically by dividing with 0.0093
    spacingsToEccDegs = rgcRFspacingsDegs/0.0093;
    sigmoidalFunction = (spacingsToEccDegs.^superGaussianExponentSigmoidalVariationExponent) ./ ...
        (spacingsToEccDegs .^superGaussianExponentSigmoidalVariationExponent + ...
         superGaussianExponentSigmoidalVariatioEccDegs50^superGaussianExponentSigmoidalVariationExponent);

    % The minimal value of gaussianPower is 2. The max  value is superGaussianMaxExponent 
    gaussianPower = 2 + (superGaussianMaxExponent-2) * (1 - sigmoidalFunction);

    ecc = sqrt(sum(rgcRFpositionsDegs.^2,2));
    eccBands = 0:0.5:100;
    for iEccBand = 1:numel(eccBands)
        idx = find (ecc >= eccBands(iEccBand)-0.25 & ecc < eccBands(iEccBand)+0.25);
        if (numel(idx > 3))
            meanGaussianPower(iEccBand) = mean(gaussianPower(idx));
            stdGaussianPower(iEccBand) = std(gaussianPower(idx));
        else
            meanGaussianPower(iEccBand) = nan;
            stdGaussianPower(iEccBand) = nan;
        end
    end

    iecc = find(~isnan(meanGaussianPower));

    hFigSuperGaussianExponents = figure(2000); clf;
    set(hFigSuperGaussianExponents, 'Name', 'superGaussianExponents');
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    theSuperGaussianExponentAxes = PublicationReadyPlotLib.generatePanelAxes(hFigSuperGaussianExponents,ff);

    scatter(theSuperGaussianExponentAxes{1,1}, ecc, gaussianPower, 32, ...
        'MarkerFaceAlpha', 0.01, ...
        'MarkerFaceColor', [0.1 0.1 0.1], ...
        'MarkerEdgeColor', 'none'); 
    hold(theSuperGaussianExponentAxes{1,1}, 'on');;
    plot(theSuperGaussianExponentAxes{1,1}, eccBands(iecc), meanGaussianPower(iecc), 'r-', 'LineWidth', 1.5);
    xlabel(theSuperGaussianExponentAxes{1,1}, 'radial eccentricity (degs)')
    ylabel(theSuperGaussianExponentAxes{1,1}, 'super Gaussian exponent');

    XLims = [0 40];
    YLims = [0 12];
    set(theSuperGaussianExponentAxes{1,1}, 'XLim', XLims, 'YLim', YLims, 'XTick', 0:5:50, 'YTick', 0:2:20);

    ff.box = 'off';
    PublicationReadyPlotLib.applyFormat(theSuperGaussianExponentAxes{1,1}, ff);
    PublicationReadyPlotLib.offsetAxes(theSuperGaussianExponentAxes{1,1}, ff, XLims, YLims, ...
        'keepXaxisSymmetric', true, ...
        'keepYaxisSymmetric', true);

    pos = get(hFigSuperGaussianExponents, 'Position');
    pos(1:2) = [600 600];
    set(hFigSuperGaussianExponents, 'Position',pos);
    drawnow;

    % Export figure
    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    figureFormat = 'png';
    %figureFormat = 'pdf';
    thePDFfileName = fullfile(theRawFiguresDir, sprintf('SuperGaussianExponents.%s', figureFormat));
    if (strcmp(figureFormat,'png'))
        NicePlot.exportFigToPNG(thePDFfileName, hFigSuperGaussianExponents,  300);
    else
        NicePlot.exportFigToPDF(thePDFfileName, hFigSuperGaussianExponents,  300);
    end

end

function [overlappingRFcenterWeightsInterpolant, xSupport, ySupport, theEllipseRotation] = generateOverlappingWeightsInterpolant(...
            allConePositions, allConeSpacings, allConeTypes, connectivityVector, ...
            exclusiveConeIndices, ...
            sensitivityAtPointOfOverlap, ...
            maxNumberOfConesOutsideContour, ...
            overlappingWeightsDivergencePattern, ...
            coneTypesIncluded, ...
            gaussianPower, ...
            theRGCindex)

    theConePositions = allConePositions(exclusiveConeIndices,:);
    theConeSpacings = allConeSpacings(exclusiveConeIndices);
    theConePoolingWeights = connectivityVector(exclusiveConeIndices);
    exclusiveConesNum = numel(theConePoolingWeights);

    spatialSupportSamples = 128;
    xSep = max(theConeSpacings)*numel(theConeSpacings);
    xx = theConePositions(:,1);
    xSupport = linspace(min(xx)-10*xSep,max(xx)+10*xSep,spatialSupportSamples);
    yy = theConePositions(:,2);
    ySupport = linspace(min(yy)-10*xSep,max(yy)+10*xSep,spatialSupportSamples);

    centerSubregionContourSamples = 32;

    if (size(theConePositions,1) >= 5)
        ellipseContourAngles = 0:20:340;
        pLevel = [];
        theSubregionContourData = mRGCMosaic.subregionEllipseFromPooledConePositions(...
                         theConePositions, centerSubregionContourSamples, ellipseContourAngles, pLevel, maxNumberOfConesOutsideContour);
        [theCenter, xAlpha, yAlpha, theRotationRadians] = RGCMosaicConstructor.helper.fit.ellipseToXYpoints( ...
            theSubregionContourData.vertices', ...
            'maxIterations', 200, ...
            'constraint', 'trace', ... % choose between {'bookstein', 'trace'}
            'nonLinear', false);
    else
        % Less than 5 cones, use subregionEllipseFromPooledCones

        theSubregionContourData = mRGCMosaic.subregionEllipseFromPooledCones(...
               	theConePositions, theConeSpacings, theConePoolingWeights, ...
                xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);

        [theCenter, xAlpha, yAlpha, theRotationRadians] = RGCMosaicConstructor.helper.fit.ellipseToXYpoints( ...
            theSubregionContourData.vertices', ...
            'maxIterations', 200, ...
            'constraint', 'trace', ... % choose between {'bookstein', 'trace'}
            'nonLinear', false);
    end


    % Overwrite center to the centroid of the exclusive cones, not that of the fitted ellipse
    theCenter = (mean(theConePositions,1))';

    % Rotation matrix
    Q = [...
        cos(theRotationRadians) -sin(theRotationRadians); ...
        sin(theRotationRadians)  cos(theRotationRadians) ...
    ];

    theEllipseRotation = theRotationRadians/pi*180;
    % Max ellipse axis
    maxAlpha = max([xAlpha yAlpha]);

    % Spatial support +/- 3xmaxAlpha
    nSamples2DGrid = 256;
    xx = linspace(-maxAlpha*5, maxAlpha*5, nSamples2DGrid);

    % Generate Gaussian weights sampler
    [Xo,Yo] = ndgrid(xx,xx);
    XY = bsxfun(@plus, [Xo(:) Yo(:)] * Q, theCenter');
    X = XY(:,1); Y = XY(:,2); 
    switch (overlappingWeightsDivergencePattern)
        case 'isotropic'
            meanAlpha = mean([xAlpha yAlpha]);
            argumentX = (X-theCenter(1))/meanAlpha;
            argumentY = (Y-theCenter(2))/meanAlpha;
        case 'inLineWithExclusiveConnections'
            argumentX = (X-theCenter(1))/xAlpha;
            argumentY = (Y-theCenter(2))/yAlpha;
        case 'orthogonalToExclusiveConnections'
            argumentX = (X-theCenter(1))/yAlpha;
            argumentY = (Y-theCenter(2))/xAlpha;
    end

    R = sqrt(argumentX.^2 + argumentY.^2);
    GaussianWeightsSampler = exp(-0.5 * R .^ gaussianPower);

    % Generate interpolant
    interpolationMethod = 'linear';
    extrapolationMethod = 'none';
    overlappingRFcenterWeightsInterpolant = ...
        griddedInterpolant(Xo+theCenter(1), Yo+theCenter(2), ...
        reshape(GaussianWeightsSampler, nSamples2DGrid*[1 1]), interpolationMethod, extrapolationMethod);


    allConeWeights = overlappingRFcenterWeightsInterpolant(allConePositions(:,1), allConePositions(:,2));
    numberOfConesWithWeightsGreaterThanSensitivityAtPointOfOverlap = numel(find(...
            (allConeWeights >= sensitivityAtPointOfOverlap) & ...
            (ismember(allConeTypes, coneTypesIncluded))...
            ));

    scalingFactor = 1;
    deltaScale = 0.01;

    Roriginal = R;
    while (numberOfConesWithWeightsGreaterThanSensitivityAtPointOfOverlap > exclusiveConesNum)
        scalingFactor = scalingFactor+deltaScale;
        % reduce the alpha and recompute interpolant
        R = Roriginal * scalingFactor;
        GaussianWeightsSampler = exp(-0.5 * R .^ gaussianPower);

        % Generate interpolant
        interpolationMethod = 'linear';
        extrapolationMethod = 'none';
        overlappingRFcenterWeightsInterpolant = ...
            griddedInterpolant(Xo+theCenter(1), Yo+theCenter(2), ...
            reshape(GaussianWeightsSampler, nSamples2DGrid*[1 1]), interpolationMethod, extrapolationMethod);

        allConeWeights = overlappingRFcenterWeightsInterpolant(allConePositions(:,1), allConePositions(:,2));
        numberOfConesWithWeightsGreaterThanSensitivityAtPointOfOverlap = numel(find(...
            (allConeWeights >= sensitivityAtPointOfOverlap) & ...
            (ismember(allConeTypes, coneTypesIncluded))...
            ));
    end

    while (numberOfConesWithWeightsGreaterThanSensitivityAtPointOfOverlap < exclusiveConesNum)
        scalingFactor = scalingFactor-deltaScale;
        % reduce the alpha and recompute interpolant
        R = Roriginal * scalingFactor;
        GaussianWeightsSampler = exp(-0.5 * R .^ gaussianPower);

        % Generate interpolant
        interpolationMethod = 'linear';
        extrapolationMethod = 'none';
        overlappingRFcenterWeightsInterpolant = ...
            griddedInterpolant(Xo+theCenter(1), Yo+theCenter(2), ...
            reshape(GaussianWeightsSampler, nSamples2DGrid*[1 1]), interpolationMethod, extrapolationMethod);

        allConeWeights = overlappingRFcenterWeightsInterpolant(allConePositions(:,1), allConePositions(:,2));
        numberOfConesWithWeightsGreaterThanSensitivityAtPointOfOverlap = numel(find(...
            (allConeWeights >= sensitivityAtPointOfOverlap) & ...
            (ismember(allConeTypes, coneTypesIncluded))...
            ));
    end

end