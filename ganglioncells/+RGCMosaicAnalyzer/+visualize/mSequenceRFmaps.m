function [sizePixels, sigmaPixels] = mSequenceRFmaps(theMRGCMosaicResponsesFullFileName, ...
        surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        sizePixels, sigmaPixels, ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsForRFmapping, ...
        residualWithRespectToNativeOpticsDefocusDiopters, ...
        generateVisualRFandConePoolingMapComboPlots, ...
        generateNearestNeighborOverlapPlots, ...
        onlyPlotPreviouslyComputedRoRincRatios, ...
        employTransparentBackground)


	if (strcmp(opticsForRFmapping, 'refractionResidualWithRespectToNativeOptics'))
		opticsInfo = sprintf('%s_residualRefraction%2.2fD',opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters);
	else
		opticsInfo = opticsForRFmapping;

	end

    if (generateNearestNeighborOverlapPlots)
        visualizeNeihboringRFmaps(theMRGCMosaicResponsesFullFileName, ...
                surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
                profileGain, zLevelsNegative, zLevelsPositive, ...
                coneFundamentalAndSpatialResolutionString, ...
                chromaticityForRFmapping, opticsInfo);
    end

	if (iscell(theMRGCMosaicResponsesFullFileName))
		[sizePixels, sigmaPixels] = visualizeMultiChromaticityMaps(theMRGCMosaicResponsesFullFileName, ...
	        surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
	        profileGain,  zLevelsNegative, zLevelsPositive, ...
            coneFundamentalAndSpatialResolutionString, ...
            chromaticityForRFmapping, opticsInfo, onlyPlotPreviouslyComputedRoRincRatios);
	else
		[sizePixels, sigmaPixels] = visualizeSingleChromaticityMaps(theMRGCMosaicResponsesFullFileName, ...
	        surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
	        profileGain, zLevelsNegative, zLevelsPositive, ...
            coneFundamentalAndSpatialResolutionString, ...
            chromaticityForRFmapping, opticsInfo, ...
            generateVisualRFandConePoolingMapComboPlots, ...
            employTransparentBackground);
	end
end


function [sizePixels, sigmaPixels] = visualizeMultiChromaticityMaps(aggregatedMRGCMosaicRFmappingResponsesFullFileNames, ...
            surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
            profileGain, zLevelsNegative, zLevelsPositive,  ...
            coneFundamentalAndSpatialResolutionString, ...
            aggregatedChromaticities, opticsInfo, onlyPlotPreviouslyComputedRoRincRatios)
 
 
    % Generate populationRoRincRatiosMatFile
    matFileName = strrep(aggregatedMRGCMosaicRFmappingResponsesFullFileNames{1}, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
    matFileName = strrep(matFileName, '.mat', 'RF.mat');
    matFileName = strrep(matFileName, '_mRGCMosaic', '');
    populationRoRincRatiosMatFile = strrep(matFileName, '_Achromatic', '');
    populationRoRincRatiosMatFile = strrep(populationRoRincRatiosMatFile, '.mat', '_RoRinc.mat');

    if (onlyPlotPreviouslyComputedRoRincRatios)
        % Save the population Ro/Rinc ratios
        load(populationRoRincRatiosMatFile, ...
            'achromaticRoRincRatios', 'centerConeIsolatingRoRincRatios', ...
            'centerConeDominances', 'centerConeNumerosities');

        targetString = 'intermediateFiles/SLIM/demos/ReidShapleyAnalyses/';
        idx = strfind(populationRoRincRatiosMatFile,targetString);
        thePDFfilename = populationRoRincRatiosMatFile(idx+numel(targetString):end);
        thePDFfilename = strrep(thePDFfilename, '.mat', '.pdf');

        % Generate RoRinc plot
        figNo = 22;
        RGCMosaicAnalyzer.visualize.centerConeIsolatingRoRincScatterPlot(figNo, ...
            achromaticRoRincRatios, centerConeIsolatingRoRincRatios, centerConeDominances, thePDFfilename);
    
        sizePixels = [];
        sigmaPixels = [];
        return;
    end

    visualizeSmoothingKernel = true;
    superimposePixelGrid = true;

    theAchromaticRFmaps = {};
    theLconeIsolatingRFmaps = {};
    theMconeIsolatingRFmaps = {};

    visualizedSensitivityRange = 0.6*[-1 1];
    [cReidShapleyLUToriginal, cReidShapleyLUTreversed] = generateReidShapleyRFmapLUT();

    % Retrieve center cone dominances
    [theRFCenterConeNumerosities, theRFCenterConeDominances] = theMRGCMosaic.allRFcenterConnectivityStats;

    % No smoothing
    sizePixels = 0;
    sigmaPixels = 0;

    % Determine smoothing kernel size based on RFmapping pixel size
    sizePixels = [];
    sigmaPixels = [];

    if (contains(coneFundamentalAndSpatialResolutionString, 'posConeFund'))
        coneFundamentalAjustmentString = 'mosaic ecc adjusted fundamentals';
    else
        coneFundamentalAjustmentString = 'foveal fundamentals'
    end

    for iChromaticity = 1:numel(aggregatedMRGCMosaicRFmappingResponsesFullFileNames)
        matFileName = strrep(aggregatedMRGCMosaicRFmappingResponsesFullFileNames{iChromaticity}, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
        matFileName = strrep(matFileName, '.mat', 'RF.mat');
        matFileName = strrep(matFileName, '_mRGCMosaic', '');

        fprintf('\nLoading computed mRGCRF RF maps and their Gaussian ellipsoid fits from %s ...', matFileName);
        load(matFileName, ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'RFmappingParamsStruct', ...
            'theRFmaps', ...
            'spatialSupportDegs');

        if (iChromaticity == 1)
            rfPixelSizeDegs = max(RFmappingParamsStruct.sizeDegs(:))/RFmappingParamsStruct.rfPixelsAcross;
            rfPixelRetinalPixelsWithin = max([1 floor(rfPixelSizeDegs/RFmappingParamsStruct.optimalResolutionDegs)]);

            sigmaPixels = 0.5*rfPixelRetinalPixelsWithin/2.0;
            sizePixels = ceil(6*sigmaPixels);
            sigmaPixelsDegs = sigmaPixels * RFmappingParamsStruct.resolutionDegs;
            sizePixelsDegs = sizePixels * RFmappingParamsStruct.resolutionDegs;
            
            if (superimposePixelGrid)
                theRFmapPixelGrid = zeros(numel(spatialSupportDegs), numel(spatialSupportDegs));
                for iRow = 1:RFmappingParamsStruct.rfPixelsAcross
                    dy = iRow*rfPixelRetinalPixelsWithin;
                    theRFmapPixelGrid(dy,:) = 1;
                end
                for iCol = 1:RFmappingParamsStruct.rfPixelsAcross
                    dx = iCol*rfPixelRetinalPixelsWithin;
                    theRFmapPixelGrid(:,dx) = 1;
                end
            else
                theRFmapPixelGrid = [];
            end
        end

        switch (aggregatedChromaticities{iChromaticity})
            case 'Achromatic'
                theAchromaticRFmaps = theRFmaps;
                theAchromaticRFmapsSpatialSupportDegs = spatialSupportDegs;
            case 'LconeIsolating'
                theLconeIsolatingRFmaps = theRFmaps;
                theLconeIsolatingRFmapsSpatialSupportDegs = spatialSupportDegs;
            case 'MconeIsolating'
                theMconeIsolatingRFmaps = theRFmaps;
                theMconeIsolatingRFmapsSpatialSupportDegs = spatialSupportDegs;
        end % switch
    end % for iChromaticity

    assert(size(theAchromaticRFmapsSpatialSupportDegs,2) == size(theLconeIsolatingRFmapsSpatialSupportDegs,2), ...
        'mismatch in RF support between achromatic and L-cone isolating RF map');

    assert(size(theAchromaticRFmapsSpatialSupportDegs,2) == size(theMconeIsolatingRFmapsSpatialSupportDegs,2), ...
        'mismatch in RF support between achromatic and M-cone isolating RF map');

    assert(theAchromaticRFmapsSpatialSupportDegs(2)-theAchromaticRFmapsSpatialSupportDegs(1) == theLconeIsolatingRFmapsSpatialSupportDegs(2)-theLconeIsolatingRFmapsSpatialSupportDegs(1), ...
         'mismatch in spatial resolution between achromatic and L-cone isolating RF map');
    assert(theAchromaticRFmapsSpatialSupportDegs(2)-theAchromaticRFmapsSpatialSupportDegs(1) == theMconeIsolatingRFmapsSpatialSupportDegs(2)-theMconeIsolatingRFmapsSpatialSupportDegs(1), ...
         'mismatch in spatial resolution between achromatic and L-cone isolating RF map');


    visualizeSmoothingKernel = true;
    if ((sizePixels == 0)||(sigmaPixels == 0))
        visualizeSmoothingKernel = false;
        fprintf('Not smoothing the measured RF\n');
        smoothingKernel = [];
    else
        fprintf('with size (degs): %f\n', sizePixelsDegs);
        fprintf('and sigma (degs): %f\n', sigmaPixelsDegs);
        smoothingKernel = fspecial('gaussian', sizePixels, sigmaPixels);
    end



    hFig = figure(111); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [0 0 2050 700]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 1, ...
        'colsNum', 5, ...
        'heightMargin',  0.01, ...
        'widthMargin',    0.02, ...
        'leftMargin',     0.05, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.04, ...
        'topMargin',      0.02);

    axConePoolingMapZoomedIn = subplot('Position', subplotPosVectors(1,1).v);
    alteredPosition =  subplotPosVectors(1,2).v;
    alteredPosition(1) = alteredPosition(1) + 0.02;
    axConePoolingMap = subplot('Position', alteredPosition);

    alteredPosition =  subplotPosVectors(1,3).v;
    alteredPosition(1) = alteredPosition(1) + 0.02;
    axAchromaticMap = subplot('Position', alteredPosition);

    alteredPosition =  subplotPosVectors(1,4).v;
    alteredPosition(1) = alteredPosition(1) + 0.01;

    axLconeIsolatingMap = subplot('Position', alteredPosition);

    alteredPosition =  subplotPosVectors(1,5).v;
    alteredPosition(1) = alteredPosition(1);
    axMconeIsolatingMap = subplot('Position', alteredPosition);

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    ffRFmap = ff;
    ffRFmap.grid = 'off';
    ffRFmap.backgroundColor = [0 0 0];

    achromaticRoRincRatios = zeros(1, numel(theTargetVisualizedRGCindices));
    centerConeIsolatingRoRincRatios = zeros(1, numel(theTargetVisualizedRGCindices));
    centerConeDominances = zeros(1, numel(theTargetVisualizedRGCindices));
    centerConeNumerosities = zeros(1, numel(theTargetVisualizedRGCindices));

    for idx = 1:numel(theTargetVisualizedRGCindices)
        fprintf('Plotting RF %d of %d\n', idx, numel(theTargetVisualizedRGCindices));
        theTargetVisualizedRGCindex = theTargetVisualizedRGCindices(idx);
        theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);

        theRFcenterConeDominance = theRFCenterConeDominances(theTargetVisualizedRGCindex);
        theRFcenterConeNumerosity = theRFCenterConeNumerosities(theTargetVisualizedRGCindex);
        theAchromaticRFmap = theAchromaticRFmaps{theTargetVisualizedRGCindex};
        theLconeIsolatingRFmap = theLconeIsolatingRFmaps{theTargetVisualizedRGCindex};
        theMconeIsolatingRFmap = theMconeIsolatingRFmaps{theTargetVisualizedRGCindex};

        % Determine plotting limits
        [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
         domainVisualizationLimits, domainVisualizationTicks, ...
         domainVisualizationLimitsSingleRF, domainVisualizationTicksSingleRF] = ...
                RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);

        % Generate mesh
        [X,Y] = meshgrid(theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(2));

        % Compute the Ro/Rpos ratio for the achromatic map
        % Ro = sum(theRFmap(:));
        % Rpos = sum(theRFmap(theRFmap>0));
        theAchromaticRoRincRatio = RoRincrementRatio(theAchromaticRFmap);
        

        % Compute the Ro/Rpos ratio for the center-cone isolating map
        if (theRFcenterConeDominance == cMosaic.LCONE_ID)
            theCenterConeIsolatingRoRincRatio = RoRincrementRatio(theLconeIsolatingRFmap);
        else
            theCenterConeIsolatingRoRincRatio = RoRincrementRatio(theMconeIsolatingRFmap);
        end

        % Save the population RoRincRatios
        achromaticRoRincRatios(idx) = theAchromaticRoRincRatio;
        centerConeIsolatingRoRincRatios(idx) = theCenterConeIsolatingRoRincRatio;
        centerConeDominances(idx) = theRFcenterConeDominance;
        centerConeNumerosities(idx) = theRFcenterConeNumerosity;

        if (~isempty(smoothingKernel))
            % Smooth the RFs
            theAchromaticRFmap = conv2(theAchromaticRFmap, smoothingKernel, 'same');
            theLconeIsolatingRFmap = conv2(theLconeIsolatingRFmap, smoothingKernel, 'same');
            theMconeIsolatingRFmap = conv2(theMconeIsolatingRFmap, smoothingKernel, 'same');
        end

        %sliceComputationMethod= 'radially symmetric average';
        sliceComputationMethod = 'sum along RF elongation axis';

        switch (sliceComputationMethod)
            case 'sum along RF elongation axis'
                % Return the rotated achromatic RF that has its elongation axis parallel to the y-axis
                % Also return the center and orientation needed for the above rotation
                [theAchromaticRFmap90, theRFmapCenter, theRFmapOrientation] = ...
                    verticallyAlignedElongationAxisRFmap(theAchromaticRFmap, [], []);

                % Apply the same transformation that we did in the achromatic RF, to rotate similarly the L- and the M-cone isolating RF maps
                theLconeIsolatingRFmap90 = verticallyAlignedElongationAxisRFmap(theLconeIsolatingRFmap, theRFmapCenter, theRFmapOrientation);
                theMconeIsolatingRFmap90 = verticallyAlignedElongationAxisRFmap(theMconeIsolatingRFmap, theRFmapCenter, theRFmapOrientation);

                % Compute slices by summing the rotated RF maps along the y-axis
                theAchromaticSlice = sum(theAchromaticRFmap90,1);
                theLconeIsolatingSlice = sum(theLconeIsolatingRFmap90,1);
                theMconeIsolatingSlice = sum(theMconeIsolatingRFmap90,1);

            case 'radially symmetric average'
 
                % Radially symmetric slices. What Reid and Shapley used
                [theAchromaticSlice, theLconeIsolatingSlice, theMconeIsolatingSlice] = ...
                    radialyAveragedSlices(theAchromaticRFmap, theLconeIsolatingRFmap, theMconeIsolatingRFmap, theRFcenterConeDominance);

            otherwise
                error('Unknown slice method: ''%s''.', sliceComputationMethod);
        end


        % Max slice
        maxSlice = max([ max(abs(theLconeIsolatingSlice(:))) max(abs(theMconeIsolatingSlice(:))) max(abs(theAchromaticSlice(:))) ]);
        ySliceOffset = 0.5*(domainVisualizationLimits(4)+domainVisualizationLimits(3))-0.2*(domainVisualizationLimits(4)-domainVisualizationLimits(3));
        ySliceScale = 0.65*(domainVisualizationLimits(4)-domainVisualizationLimits(3));

        % Shift so we can display the slices on the same plot as the RF maps
        theAchromaticSlice = ySliceOffset  + theAchromaticSlice / maxSlice * ySliceScale;
        theLconeIsolatingSlice = ySliceOffset  + theLconeIsolatingSlice / maxSlice * ySliceScale;
        theMconeIsolatingSlice = ySliceOffset  + theMconeIsolatingSlice / maxSlice * ySliceScale;

        if (~isempty(theRFmapPixelGrid))
            % Superimpose theRFmapPixelGrid
            theAchromaticRFmap(theRFmapPixelGrid==1) = 0;
            theLconeIsolatingRFmap(theRFmapPixelGrid==1) = 0;
            theMconeIsolatingRFmap(theRFmapPixelGrid==1) = 0;
        end
        
        maxAllMaps = max([max(abs(theAchromaticRFmap(:))) max(abs(theLconeIsolatingRFmap(:))) max(abs(theMconeIsolatingRFmap(:)))]);
        cla(axAchromaticMap);
        cla(axLconeIsolatingMap);
        cla(axMconeIsolatingMap);
        imagesc(axAchromaticMap, ...
            theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ...
            theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(2), ...
            theAchromaticRFmap/maxAllMaps);

       
        if (visualizeSmoothingKernel) && (~isempty(smoothingKernel))
            [~,idx] = max(abs(theAchromaticRFmap(:)));
            [deltaFunctionRow, deltaFunctionCol] = ind2sub(size(theAchromaticRFmap), idx);
            theDeltaFunction = zeros(size(theAchromaticRFmap));;
            % Place the smoothing kernel up and to the right
            ddX = round(deltaFunctionCol-3*size(theAchromaticRFmap,1)/10);
            ddY = round(deltaFunctionRow+2*size(theAchromaticRFmap,1)/10);
            xo = max([2, min([size(theAchromaticRFmap,1)-2 ddX])]);
            yo = max([2 min([size(theAchromaticRFmap,1)-2 ddY])]);
            theDeltaFunction(yo, xo) = 1;
            theVisualizedSmoothingKernel = conv2(theDeltaFunction, smoothingKernel, 'same');
            hold(axAchromaticMap, 'on');
            % 1, 1 x sigma, 5%
            zLevelsPeakOneSigmaTwoSigma = [1 exp(-0.5*(1)^2)  5/100];
            contourf(axAchromaticMap, X,Y, theVisualizedSmoothingKernel/max(theVisualizedSmoothingKernel(:)),  ...
                zLevelsPeakOneSigmaTwoSigma, 'Color', 'none', 'LineWidth', 1.0, 'LineStyle', '-', 'LineColor', [1 1 1]);
        end

        % Slice line config
        outerSliceColor = [1 0.7 0.4];
        outerSliceWidth = 4;
        innerSliceColor = [0.7 0.5 0.0];
        innerSliceWidth = 1.5;

        % Add the 1D slice
        if (~isempty(smoothingKernel))
            plot(axAchromaticMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ySliceOffset + 0*theAchromaticSlice, ...
                '--', 'LineWidth', 1.0, 'Color', innerSliceColor);
            plot(axAchromaticMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theAchromaticSlice, ...
                '-', 'LineWidth', outerSliceWidth, 'Color', outerSliceColor);
            plot(axAchromaticMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theAchromaticSlice, ...
                '-', 'LineWidth', innerSliceWidth, 'Color', innerSliceColor);
        end
        set(axAchromaticMap, 'XTickLabel', [], 'YTickLabel', []);

        title(axAchromaticMap, sprintf('achromatic RF (Ro/Rinc: %2.3f)\n(%s)\n', ...
                theAchromaticRoRincRatio, ...
                coneFundamentalAjustmentString));

        imagesc(axLconeIsolatingMap, ...
            theLconeIsolatingRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ...
            theLconeIsolatingRFmapsSpatialSupportDegs+spatialSupportCenterDegs(2), ...
            theLconeIsolatingRFmap/maxAllMaps);

        % Add the 1D slice
        if (~isempty(smoothingKernel))
            hold(axLconeIsolatingMap, 'on');
            plot(axLconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ySliceOffset + 0*theAchromaticSlice, ...
                '--', 'LineWidth', 1.0, 'Color', innerSliceColor);
            plot(axLconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theLconeIsolatingSlice, ...
                '-', 'LineWidth', outerSliceWidth, 'Color', outerSliceColor);
            plot(axLconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theLconeIsolatingSlice, ...
                '-', 'LineWidth', innerSliceWidth, 'Color', innerSliceColor);
        end

        set(axLconeIsolatingMap, 'XTickLabel', [], 'YTickLabel', []);

        if (theRFcenterConeDominance == cMosaic.LCONE_ID)
            title(axLconeIsolatingMap, sprintf('L-cone isolating RF (Ro/Rinc: %2.3f)\n(%s)\n', ...
                theCenterConeIsolatingRoRincRatio, ...
                coneFundamentalAjustmentString));
        else
            title(axLconeIsolatingMap, sprintf('L-cone isolating RF\n(%s)\n',coneFundamentalAjustmentString));
        end

        imagesc(axMconeIsolatingMap, ...
            theMconeIsolatingRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ...
            theMconeIsolatingRFmapsSpatialSupportDegs+spatialSupportCenterDegs(2), ...
            theMconeIsolatingRFmap/maxAllMaps);

        % Add the 1D slice
        if (~isempty(smoothingKernel))
            hold(axMconeIsolatingMap, 'on');
            plot(axMconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), ySliceOffset + 0*theAchromaticSlice, ...
                '--', 'LineWidth', 1.0, 'Color', innerSliceColor);
            plot(axMconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theMconeIsolatingSlice, ...
                '-', 'LineWidth', outerSliceWidth, 'Color', outerSliceColor);
            plot(axMconeIsolatingMap, ...
                theAchromaticRFmapsSpatialSupportDegs+spatialSupportCenterDegs(1), theMconeIsolatingSlice, ...
                '-', 'LineWidth', innerSliceWidth, 'Color', innerSliceColor);
        end

        if (theRFcenterConeDominance == cMosaic.MCONE_ID)
            title(axMconeIsolatingMap, sprintf('M-cone isolating RF (Ro/Rinc: %2.3f)\n(%s)\n', ...
                theCenterConeIsolatingRoRincRatio, ...
                coneFundamentalAjustmentString));
        else
            title(axMconeIsolatingMap, sprintf('M-cone isolating RF\n(%s)\n',coneFundamentalAjustmentString));
        end


        set(axAchromaticMap, 'CLim', visualizedSensitivityRange);
        colormap(axAchromaticMap, cReidShapleyLUToriginal);
        axis(axAchromaticMap, 'image'); axis(axAchromaticMap, 'xy');
        RGCMosaicAnalyzer.visualize.setLimsAndTicks(...
            axAchromaticMap, domainVisualizationTicks, domainVisualizationLimits);
        set(axAchromaticMap, 'XTickLabel', [], 'YTickLabel', []);

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(axAchromaticMap,ffRFmap);

        set(axLconeIsolatingMap, 'CLim', visualizedSensitivityRange);
        colormap(axLconeIsolatingMap, cReidShapleyLUToriginal);
        axis(axLconeIsolatingMap, 'image'); axis(axLconeIsolatingMap, 'xy');
        RGCMosaicAnalyzer.visualize.setLimsAndTicks(...
            axLconeIsolatingMap, domainVisualizationTicks, domainVisualizationLimits);
        set(axLconeIsolatingMap, 'XTickLabel', [], 'YTickLabel', []);

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(axLconeIsolatingMap,ffRFmap);

        set(axMconeIsolatingMap, 'CLim', visualizedSensitivityRange);
        colormap(axMconeIsolatingMap, cReidShapleyLUToriginal);
        axis(axMconeIsolatingMap, 'image'); axis(axMconeIsolatingMap, 'xy');
        RGCMosaicAnalyzer.visualize.setLimsAndTicks(...
            axMconeIsolatingMap, domainVisualizationTicks, domainVisualizationLimits);
        set(axMconeIsolatingMap, 'XTickLabel', [], 'YTickLabel', []);

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(axMconeIsolatingMap,ffRFmap);

        pdfFileName = sprintf('RFs_%s_%s_RGC_%d.pdf', opticsInfo, coneFundamentalAndSpatialResolutionString, theTargetVisualizedRGCindex);
        figNo = [];
        scaleBarDegs = 7.5/60;
        RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(figNo, ...
                theMRGCMosaic, theTargetVisualizedRGCindex, '', ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
                'fixedScaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'noGrid', true, ...
                'renderLineWeightingFunctionPlots', false, ...
                'plotTitle', sprintf('cone pooling map\n'), ...
                'figureHandle', hFig, ...
                'axesHandle', axConePoolingMap, ...
                'figureFormat', ff);
        ylabel(axConePoolingMap, '');
        

        RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(figNo, ...
                theMRGCMosaic, theTargetVisualizedRGCindex, pdfFileName, ...
                'domainVisualizationLimits', domainVisualizationLimitsSingleRF, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
                'fixedScaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'noGrid', true, ...
                'renderLineWeightingFunctionPlots', false, ...
                'plotTitle', sprintf('cone pooling map (zoomed-in)\n'), ...
                'figureHandle', hFig, ...
                'axesHandle', axConePoolingMapZoomedIn, ...
                'figureFormat', ff);
    end % idx  

    % Save the population Ro/Rinc ratios
    save(populationRoRincRatiosMatFile, ...
            'achromaticRoRincRatios', 'centerConeIsolatingRoRincRatios', ...
            'centerConeDominances', 'centerConeNumerosities');

    % Generate RoRinc plot
    figNo = 22;
    RGCMosaicAnalyzer.visualize.centerConeIsolatingRoRincScatterPlot(figNo, ...
        achromaticRoRincRatios, centerConeIsolatingRoRincRatios, centerConeDominances, '');
    
end


function visualizeNeihboringRFmaps(theMRGCMosaicResponsesFullFileName, ...
        surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive,  ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsInfo)
    matFileName = strrep(theMRGCMosaicResponsesFullFileName, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
    matFileName = strrep(matFileName, '.mat', 'RF.mat');
    matFileName = strrep(matFileName, '_mRGCMosaic', '');

    fprintf('\nLoading computed mRGCRF RF maps and their Gaussian ellipsoid fits from %s ...', matFileName);
    load(matFileName, ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'RFmappingParamsStruct', ...
            'theRFmaps', ...
            'spatialSupportDegs', ...
            'theFittedGaussianEllipsoids');

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';

    %sliceComputationMethod = 'radially symmetric average';
    sliceComputationMethod = 'sum along RF elongation axis';

    for idx = 1:numel(theTargetVisualizedRGCindices)
        % Retrieve the RGC index 
        theTargetVisualizedRGCindex = theTargetVisualizedRGCindices(idx);

        % Retrieve the RF map
        theTargetRFmap = theRFmaps{theTargetVisualizedRGCindex};
        theTargetRFmap = theTargetRFmap / max(theTargetRFmap(:));

        % Retrieve the Gaussian ellipsoid fit to the positive part of the RFmap
        theFittedEllipsoid  = theFittedGaussianEllipsoids{theTargetVisualizedRGCindex};
        theFittedTargetRFmap = theFittedEllipsoid.rfMap;

        % Compute profile
        switch (sliceComputationMethod)
            case 'sum along RF elongation axis'
                % Return the rotated achromatic RF that has its elongation axis parallel to the y-axis
                % Also return the center and orientation needed for the above rotation
                theTargetRFmap90 = verticallyAlignedElongationAxisRFmap(theTargetRFmap, [], []);
                % Compute slices by summing the rotated RF maps along the y-axis
                theTargetRFprofile = sum(theTargetRFmap90,1);
                % Compute slices by summing the rotated RF maps along the x-axis
                theTargetRFprofile2 = sum(theTargetRFmap90,2);

            case 'radially symmetric average'
                % Radially symmetric slices.
                theTargetRFprofile = circularlySymmetricSlice(theTargetRFmap, []);

            otherwise
                error('Unknown slice method: ''%s''.', sliceComputationMethod);
        end

        % Find the closest RGC index
        theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);
        [~, idx] = MosaicConnector.pdist2(theMRGCMosaic.rgcRFpositionsDegs, ...
            theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:), ...
            'smallest', 2);
        theNearestVisualizedRGCindex = idx(2);

        % Retrieve the RF map
        theNearestRGCRFmap = theRFmaps{theNearestVisualizedRGCindex};
        theNearestRGCRFmap = theNearestRGCRFmap / max(theNearestRGCRFmap(:));

        % Retrieve the Gaussian ellipsoid fit to the positive part of the RFmap
        theFittedEllipsoid  = theFittedGaussianEllipsoids{theNearestVisualizedRGCindex};
        theFittedNearestRFmap = theFittedEllipsoid.rfMap;

        % Compute profile
        switch (sliceComputationMethod)
            case 'sum along RF elongation axis'
                % Return the rotated achromatic RF that has its elongation axis parallel to the y-axis
                % Also return the center and orientation needed for the above rotation
                theTargetRFmap90 = verticallyAlignedElongationAxisRFmap(theNearestRGCRFmap, [], []);
                % Compute slices by summing the rotated RF maps along the y-axis
                theNearestRFprofile = sum(theTargetRFmap90,1);
                % Compute slices by summing the rotated RF maps along the y-axis
                theNearestRFprofile2 = sum(theTargetRFmap90,2);
                maxAll2  = max([max(theTargetRFprofile2) max(theNearestRFprofile2)]);

            case 'radially symmetric average'
                % Radially symmetric slices.
                theNearestRFprofile = circularlySymmetricSlice(theNearestRGCRFmap, []);
                maxAll2 = 0;

            otherwise
                error('Unknown slice method: ''%s''.', sliceComputationMethod);
        end

        maxAllY = max([max(theTargetRFprofile) max(theNearestRFprofile)]);
        maxAll = max([maxAllY maxAll2]);

        ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');


        % Plot the target RF map, along with its profile, and the RFmap contour
        pdfFileName = sprintf('%s_%stargetRFmap%d.pdf', opticsInfo, chromaticityForRFmapping, theTargetVisualizedRGCindex);
            
        figNo = theTargetVisualizedRGCindex+31000;
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        % Visualize the target RF map, along with the contour at 1 sigma of the target RF and of the nearest RGC RF
        imagesc(ax, spatialSupportDegs, spatialSupportDegs, theTargetRFmap);
        axis(ax, 'image'); axis(ax, 'xy')
        set(ax, 'CLim', max(theTargetRFmap(:))*[-1 1]);
        hold (ax, 'on');
        contour(ax, theFittedEllipsoid.xSupportDegs, theFittedEllipsoid.ySupportDegs, theFittedTargetRFmap/max(theFittedTargetRFmap(:)), [exp(-0.502) exp(-0.50)], ...
                'Color', [1 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');
        contour(ax, theFittedEllipsoid.xSupportDegs, theFittedEllipsoid.ySupportDegs, theFittedNearestRFmap/max(theFittedNearestRFmap(:)), [exp(-0.502) exp(-0.50)], ...
                'Color', [0 0 1], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');

        % Plot the profile
        plot(ax, spatialSupportDegs, ...
            theTargetRFprofile/maxAllY*0.5*(spatialSupportDegs(end)-spatialSupportDegs(1))+theFittedEllipsoid.y0, ... 
            'r-', 'LineWidth', 1.5);
        hold(ax, 'off');
        colormap(ax, gray)

        xlabel(ax, 'position, x (degs)');
        ylabel(ax, 'positions, y (degs)');

        ff.box = 'off';
        ff.tickDir = 'in';
        ff.grid = 'off';
        ff.backgroundColor = [1 1 1];

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        % Export the PDF
        thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
        % Using the export_fig which supports transparent background
        export_fig(hFig, thePDFfileName, '-pdf', '-native');


        % Plot the nearest to the target RF map, along with its profile, and the RFmap contour
        pdfFileName = sprintf('%s_%snearestRFmap%d.pdf', opticsInfo, chromaticityForRFmapping, theTargetVisualizedRGCindex);
            
        figNo = theTargetVisualizedRGCindex+32000;
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        imagesc(ax, spatialSupportDegs, spatialSupportDegs, theNearestRGCRFmap);
        axis(ax, 'image'); axis(ax, 'xy')
        set(ax, 'CLim', max(theNearestRGCRFmap(:))*[-1 1]);
        hold (ax, 'on');
        contour(ax, theFittedEllipsoid.xSupportDegs, theFittedEllipsoid.ySupportDegs, ...
            theFittedTargetRFmap/max(theFittedTargetRFmap(:)), [exp(-0.502) exp(-0.50)], ...
                'Color', [1 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');

        contour(ax, theFittedEllipsoid.xSupportDegs, theFittedEllipsoid.ySupportDegs, ...
            theFittedNearestRFmap/max(theFittedNearestRFmap(:)), [exp(-0.502) exp(-0.50)], ...
                'Color', [0 0 1], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');

        plot(ax, spatialSupportDegs, theNearestRFprofile/maxAllY*0.5*(spatialSupportDegs(end)-spatialSupportDegs(1))+theFittedEllipsoid.y0, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');
        colormap(ax, gray)

        xlabel(ax, 'position, x (degs)');
        ylabel(ax, 'positions, y (degs)');

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        % Export the PDF
        thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
        % Using the export_fig which supports transparent background
        export_fig(hFig, thePDFfileName, '-pdf', '-native');


        % Plot the target and nearest RF map x-profiles 
        pdfFileName = sprintf('%s_%soverlapsXRFmap%d.pdf', opticsInfo, chromaticityForRFmapping, theTargetVisualizedRGCindex);
            
        figNo = theTargetVisualizedRGCindex+33000;
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        plot(ax, spatialSupportDegs, theTargetRFprofile/maxAll, 'r-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, spatialSupportDegs, theNearestRFprofile/maxAll, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');

        set(ax, 'YLim', [-1 1]);
        xlabel(ax, 'position, x (degs)');
        ylabel(ax, 'sensitivity');

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        % Export the PDF
        thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
        % Using the export_fig which supports transparent background
        export_fig(hFig, thePDFfileName, '-pdf', '-native');



        % Plot the target and nearest RF map y-profiles 
        pdfFileName = sprintf('%s_%soverlapsYRFmap%d.pdf', opticsInfo, chromaticityForRFmapping, theTargetVisualizedRGCindex);
            
        figNo = theTargetVisualizedRGCindex+34000;
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        plot(ax, spatialSupportDegs, theTargetRFprofile2/maxAll, 'r-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, spatialSupportDegs, theNearestRFprofile2/maxAll, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');

        set(ax, 'YLim', [-1 1]);
        xlabel(ax, 'position, y (degs)');
        ylabel(ax, 'sensitivity');

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        % Export the PDF
        thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
        % Using the export_fig which supports transparent background
        export_fig(hFig, thePDFfileName, '-pdf', '-native');

        sprintf('Hit enter to continue with cell %d of %d\n', idx, numel(theTargetVisualizedRGCindices));
        pause
    end % idx

end


function [sizePixels, sigmaPixels] = visualizeSingleChromaticityMaps(theMRGCMosaicResponsesFullFileName, ...
        surroundConnectedParamsStruct, theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive,  ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsInfo, ...
        generateVisualRFandConePoolingMapComboPlots, ...
        employTransparentBackground)

	matFileName = strrep(theMRGCMosaicResponsesFullFileName, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
    matFileName = strrep(matFileName, '.mat', 'RF.mat');
    matFileName = strrep(matFileName, '_mRGCMosaic', '');

    fprintf('\nLoading computed mRGCRF RF maps and their Gaussian ellipsoid fits from %s ...', matFileName);
    load(matFileName, ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'RFmappingParamsStruct', ...
            'theRFmaps', ...
            'spatialSupportDegs', ...
            'theFittedGaussianEllipsoids');

    visualizeSmoothingKernel = false;
    switch (chromaticityForRFmapping)
        case 'Achromatic'
            theProfileColor = RGCMosaicConstructor.constants.achromaticColor;
            visualizeSmoothingKernel = true;
        case 'LconeIsolating'
            theProfileColor = RGCMosaicConstructor.constants.LcenterColor;
        case 'MconeIsolating'
            theProfileColor = RGCMosaicConstructor.constants.McenterColor;
    end % switch

    radialEccDegs = zeros(1, numel(theTargetVisualizedRGCindices));
    radialEccMMs = zeros(1, numel(theTargetVisualizedRGCindices));
    majorSigmaDegs = zeros(1, numel(theTargetVisualizedRGCindices));
    minorSigmaDegs = zeros(1, numel(theTargetVisualizedRGCindices));
    sigmaDegs= zeros(1, numel(theTargetVisualizedRGCindices));
    majorSigmaMicrons = zeros(1, numel(theTargetVisualizedRGCindices));
    minorSigmaMicrons = zeros(1, numel(theTargetVisualizedRGCindices));
    sigmaMicrons = zeros(1, numel(theTargetVisualizedRGCindices));

    % All the fitted Gaussian ellipsoids
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';

    % Plot the fitted Gaussian ellipsoids
    pdfFileName = sprintf('%s_%s_fittedGaussianEllipsoids.pdf', opticsInfo, chromaticityForRFmapping);

    figNo = 200;
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    if (employTransparentBackground)
        set(hFig, 'Color', 'none');
    end 
    ax = theAxes{1,1};
    hold(ax, 'on');

    for idx = 1:numel(theTargetVisualizedRGCindices)
        fprintf('Plotting RF %d of %d\n', idx, numel(theTargetVisualizedRGCindices));
        theTargetVisualizedRGCindex = theTargetVisualizedRGCindices(idx);
        theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);
        theRFmap = theRFmaps{theTargetVisualizedRGCindex};
        theFittedEllipsoid  = theFittedGaussianEllipsoids{theTargetVisualizedRGCindex};
        if (isempty(theFittedEllipsoid))
            fprintf(2,'*** No fit for cell %d\n', theTargetVisualizedRGCindex);
            continue
        end

        theFittedEllipsoidMap = theFittedEllipsoid.rfMap;
        theFittedEllipsoidMap = theFittedEllipsoidMap / max(theFittedEllipsoidMap(:));

        radialEccDegs(idx) = sqrt(sum(theRGCpositionDegs.^2,2));
        radialEccMMs(idx) = radialEccDegs(idx) * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;
        majorSigmaDegs(idx) = theFittedEllipsoid.majorSigmaDegs;
        minorSigmaDegs(idx) = theFittedEllipsoid.minorSigmaDegs;
        sigmaDegs(idx) = sqrt(majorSigmaDegs(idx)*minorSigmaDegs(idx));
        majorSigmaMicrons(idx) = theFittedEllipsoid.majorSigmaMicrons;
        minorSigmaMicrons(idx) = theFittedEllipsoid.minorSigmaMicrons;
        sigmaMicrons(idx) = sqrt(majorSigmaMicrons(idx)*minorSigmaMicrons(idx));

        supportDegsX = theFittedEllipsoid.xSupportDegs+theMRGCMosaic.eccentricityDegs(1);
        supportMMsX = supportDegsX  * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;

        supportDegsY = theFittedEllipsoid.ySupportDegs+theMRGCMosaic.eccentricityDegs(2);
        supportMMsY = supportDegsY  * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;

        contourf(ax, supportMMsX, supportMMsY, ...
                theFittedEllipsoidMap, [exp(-0.502) exp(-0.50)], ...
                'Color', [0 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');
        drawnow;
    end

    XLimsMMs = round((theMRGCMosaic.eccentricityDegs(1) + 0.6*theMRGCMosaic.sizeDegs(1)*[-1 1]) * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3 * 10)/10;
    YLimsMMs = round((theMRGCMosaic.eccentricityDegs(2) + 0.6*theMRGCMosaic.sizeDegs(2)*[-1 1])* theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3 * 10)/10;

    if (XLimsMMs(1) == XLimsMMs(2))
        XLimsMMs = XLimsMMs(1) + 0.6*theMRGCMosaic.sizeDegs(1)*[-1 1] * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;
    end
    if (YLimsMMs(1) == YLimsMMs(2))
        YLimsMMs = YLimsMMs(1) + 0.6*theMRGCMosaic.sizeDegs(2)*[-1 1] * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;
    end

    xTicksMMs = round(theMRGCMosaic.eccentricityDegs(1) * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3 * 10)/10 + -1:0.2:1;
    yTicksMMs = round(theMRGCMosaic.eccentricityDegs(2) * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3 * 10)/10 + -1:0.2:1;
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'XLim', XLimsMMs, 'YLim', YLimsMMs, 'XTick', xTicksMMs, 'YTick', yTicksMMs);

    % 25 microns size reference
    plot(ax, XLimsMMs(1) + (20+[0 25])/1000, YLimsMMs(1) + 30/1000*[1 1], 'k-', 'LineWidth', 1.5);

    % 500 microns size reference
    plot(ax, XLimsMMs(2) - (20+[0 500])/1000, YLimsMMs(1) + 30/1000*[1 1], 'k-', 'LineWidth', 1.5);

    xlabel(ax, 'eccentricity, x (mm)');
    ylabel(ax, 'eccentricity, y (mm)');

    ff.box = 'off';
    ff.tickDir = 'in';
    ff.grid = 'off';
    if (employTransparentBackground)
        ff.backgroundColor = 'none';
    end

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    
    % Export the PDF
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    if (employTransparentBackground)
        % Using the export_fig which supports transparent background
        export_fig(hFig, thePDFfileName, '-pdf', '-native');
    else 
        NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);
    end


    % Plot the fitted Gaussian radii as a function of eccentricity together with EJ's data
    pdfFileName = sprintf('%s_%s_RFdiametersAsAFunctionOfEccentricity.pdf', opticsInfo, chromaticityForRFmapping);


    displayMinorAndMajorAxisDiameters = false;
	employLogYaxis = false;
	employLogXaxis = ~true;
    RGCMosaicAnalyzer.visualize.mosaicRetinalDiametersAgainstMacaqueInVitroRFdiameters(201, ff, ...
        radialEccMMs, minorSigmaMicrons, majorSigmaMicrons, sigmaMicrons, ...
        displayMinorAndMajorAxisDiameters, employLogYaxis, employLogXaxis, ...
        pdfFileName);

    % Automatic smoothing
    %sizePixels = []; sigmaPixels = [];

    % No smoothing
    sizePixels = 0; sigmaPixels = 0;

    if (generateVisualRFandConePoolingMapComboPlots)
        for idx = 1:numel(theTargetVisualizedRGCindices)
            theTargetVisualizedRGCindex = theTargetVisualizedRGCindices(idx);
            theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);
            theRFmap = theRFmaps{theTargetVisualizedRGCindex};
            theFittedEllipsoid  = theFittedGaussianEllipsoids{theTargetVisualizedRGCindex};
            
            % Determine plotting limits
            [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
            domainVisualizationLimits, domainVisualizationTicks, ...
            domainVisualizationLimitsSingleRF, domainVisualizationTicksSingleRF] = ...
                RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);


            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');

            % Plot the computed RF map
            pdfFileName = sprintf('%s_%sRFmap%d.pdf', opticsInfo, chromaticityForRFmapping, theTargetVisualizedRGCindex);
            
            figNo = theTargetVisualizedRGCindex+10000;
            hFig = figure(figNo); clf;
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};
            
       		[sizePixels, sigmaPixels, spatialSupportCenterDegs] = visualizeRFMap(theMRGCMosaic, ...
                theTargetVisualizedRGCindex, spatialSupportDegs, theRFmap, theFittedEllipsoid, ...
                zLevelsNegative, zLevelsPositive, theProfileColor, profileGain, ...
                visualizeSmoothingKernel, sizePixels, sigmaPixels, ...
                hFig, ax, ff, pdfFileName, chromaticityForRFmapping);

            figNo = theTargetVisualizedRGCindex+20000;
            hFig = figure(figNo); clf;
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};

            % Plot the cone weights map
            pdfFileName = sprintf('conePoolingWeightsMap%d.pdf', theTargetVisualizedRGCindex);
            RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(figNo, ...
                theMRGCMosaic, theTargetVisualizedRGCindex, pdfFileName, ...
                'domainVisualizationLimits', domainVisualizationLimitsSingleRF, ...
                'domainVisualizationTicks', domainVisualizationTicksSingleRF, ...
                'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
                'fixedScaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'noGrid', true, ...
                'plotTitle', '', ...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'figureFormat', ff);

            % Pause and ask user to unpause
            RGCMosaicConstructor.helper.queryUserFor.unpausingExecution(sprintf('RFmap %d of %d\n', idx, numel(theTargetVisualizedRGCindices)));
        
        end % for idx
   end  %  if (generateVisualRFandConePoolingMapComboPlots)
end


function visualizeConeWeightsRFmap(theMRGCMosaic, theTargetVisualizedRGCindex, spatialSupportDegs, hFig, ax, ff, thePDFfileName)

    theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);

    % Determine plotting limits
    [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
     domainVisualizationLimits, domainVisualizationTicks, ...
     domainVisualizationLimitsSingleRF, domainVisualizationTicksSingleRF] = ...
        RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);

    RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(0, ...
            theMRGCMosaic, theTargetVisualizedRGCindex, thePDFfileName, ...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'figureFormat', ff, ...
            'noGrid', true, ...
            'domainVisualizationLimits', domainVisualizationLimitsSingleRF, ...
            'domainVisualizationTicks', domainVisualizationTicksSingleRF, ...
            'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
            'fixedScalaBarDegs', scaleBarDegs, ...
            'doNotLabelScaleBar', true, ...
            'plotTitle', '');
end


function [sizePixels, sigmaPixels, spatialSupportCenterDegs] = visualizeRFMap(theMRGCMosaic, ...
    theTargetVisualizedRGCindex, theRFmapSpatialSupportDegs, theRFmap,  theFittedEllipsoid, ...
    zLevelsNegative, zLevelsPositive, theProfileColor, profileGain, ...
    visualizeSmoothingKernel, sizePixels, sigmaPixels, ...
    hFig, ax, ff, pdfFileName, plotTitle)

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';

    theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);
    idx = theMRGCMosaic.singleCellConnectivityStats(theTargetVisualizedRGCindex, 'center', ...
            'inputConeIndicesOnly', true);
    theConeApertureDiameterDegs = mean(theMRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(idx));

    [sizePixelsAutomatic, sigmaPixelsAutomatic, sizePixelsAutomaticDegs, sigmaPixelsAutomaticDegs] = ...
        smoothingKernelFromConeAperture(theConeApertureDiameterDegs, theRFmapSpatialSupportDegs);

    [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
        domainVisualizationLimitsFull, domainVisualizationTicksFull, ...
        domainVisualizationLimits, domainVisualizationTicks] = ...
            RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);


    if (isempty(sizePixels) || (isempty(sigmaPixels)))
        % Set smoothing kernel size based on cone aperture (automatic)
        sizePixels = sizePixelsAutomatic;
        sigmaPixels = sigmaPixelsAutomatic;
        fprintf('Employing automatic RF smoothing kernel\n');
        fprintf('with size (degs): %f\n', sizePixelsAutomaticDegs);
        fprintf('and sigma (degs): %f\n', sigmaPixelsAutomaticDegs);
        % Smooth the RF for the contours
        smoothedRFmap = conv2(theRFmap, fspecial('gaussian', sizePixels, sigmaPixels), 'same');
    else
        if ((sizePixels == 0)||(sigmaPixels == 0))
            smoothedRFmap = theRFmap;
            visualizeSmoothingKernel = false;
            fprintf('Not smoothing the measured RF\n');
        else
            dx = sigmaPixelsAutomaticDegs/sigmaPixelsAutomatic;
            fprintf('Employing custom RF smoothing kernel\n');
            fprintf('with size (degs): %f\n', sizePixels * dx);
            fprintf('and sigma (degs): %f\n', sigmaPixels * dx);
            % Smooth the RF for the contours
            smoothedRFmap = conv2(theRFmap, fspecial('gaussian', sizePixels, sigmaPixels), 'same');
        end
    end

    if (visualizeSmoothingKernel)
        [~,idx] = max(abs(theRFmap(:)));
        [deltaFunctionRow, deltaFunctionCol] = ind2sub(size(theRFmap), idx);
        theDeltaFunction = theRFmap*0;
        dd = round(size(theRFmap,1)/10);
        theDeltaFunction(min([size(theRFmap,1)-2 deltaFunctionRow+dd]), min([size(theRFmap,1)-2 deltaFunctionCol+dd])) = 1;
        smoothingKernel = conv2(theDeltaFunction, fspecial('gaussian', sizePixels, sigmaPixels), 'same');
    end

    % Generate mesh
    [X,Y] = meshgrid(theRFmapSpatialSupportDegs+spatialSupportCenterDegs(1), theRFmapSpatialSupportDegs+spatialSupportCenterDegs(2));
    
    imagesc(ax, theRFmapSpatialSupportDegs+spatialSupportCenterDegs(1), theRFmapSpatialSupportDegs+spatialSupportCenterDegs(2), ...
        theRFmap/max(abs(theRFmap(:))));
    hold(ax, 'on');

    if ((sizePixels == 0)||(sigmaPixels == 0))
    else
        % The RF contours
        contour(ax, X,Y, smoothedRFmap, zLevelsNegative, 'Color', [0 0 1], 'LineWidth', 2.0, 'LineStyle', ':');
        contour(ax, X,Y, smoothedRFmap, zLevelsPositive, 'Color', [1 0 0], 'LineWidth', 1.5, 'LineStyle', '-');

        % The smoothing kernel
        if (visualizeSmoothingKernel)
            contourf(ax, X,Y, smoothingKernel/max(smoothingKernel(:)),  0.05:0.2:1.0, 'Color', 'none', 'LineWidth', 1.0, 'LineStyle', '-');
            set(ax, 'CLim', [0 1]);
        end
    end

    if (~isempty(theFittedEllipsoid))

        theFittedEllipsoidMap = theFittedEllipsoid.rfMap;
        theFittedEllipsoidMap = theFittedEllipsoidMap / max(theFittedEllipsoidMap(:));

        hold(ax, 'on');
        contour(ax,theFittedEllipsoid.xSupportDegs+spatialSupportCenterDegs(1), ...
                        theFittedEllipsoid.ySupportDegs+spatialSupportCenterDegs(2), ...
                theFittedEllipsoidMap, [exp(-0.5) exp(-0.51)], ...
                'Color', [1 1 1], 'LineWidth', 4.0, 'LineStyle', '-');
        contour(ax,theFittedEllipsoid.xSupportDegs+spatialSupportCenterDegs(1), ...
                        theFittedEllipsoid.ySupportDegs+spatialSupportCenterDegs(2), ...
                theFittedEllipsoidMap, [exp(-0.5) exp(-0.51)], ...
                'Color', [0 0 0], 'LineWidth', 2.0, 'LineStyle', '-');
    end

    amplitude = domainVisualizationLimits(4)-domainVisualizationLimits(3);
    baselineX = domainVisualizationLimits(3) + amplitude*0.1;
    baselineY = domainVisualizationLimits(1) + amplitude*0.1;

    horizontalLineWeightingFunction = 0.2*profileGain*sum(smoothedRFmap,1);
    verticalLineWeightingFunction = 0.2*profileGain*sum(smoothedRFmap,2);

    
    y = baselineX+0.5*amplitude*horizontalLineWeightingFunction;
    x = theRFmapSpatialSupportDegs+spatialSupportCenterDegs(1);
    %RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax, x, y, baselineX, theProfileColor, theProfileColor*0.5, 0.7, 0.1);
    plot(ax, x,y, 'k-', 'Color', [1 0 0], 'LineWidth', 2.0);
    plot(ax, x, baselineX + x*0, 'k--', 'Color', [0 0 0], 'LineWidth', 1.0);


    y = theRFmapSpatialSupportDegs+spatialSupportCenterDegs(2);
    x = baselineY+0.5*amplitude*verticalLineWeightingFunction;
    %RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax, x, y, baselineY, theProfileColor, theProfileColor*0.5, 0.7, 0.1);
    plot(ax, x,y, 'k-', 'Color', [1 0 0], 'LineWidth', 2.0);
    plot(ax, baselineY + y*0,y, 'k--', 'Color', [0 0 0], 'LineWidth', 1.0);

    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'CLim', 0.2*[-1 1]);
    colormap(ax, brewermap(1024, '*RdBu'));

    RGCMosaicAnalyzer.visualize.setLimitsAndTicks(...
        ax, domainVisualizationTicks, domainVisualizationLimits);

    title(ax, plotTitle);

    ff.box = 'on';
    ff.tickDir = 'in';
    ff.grid = 'off';

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end


function [sizePixels, sigmaPixels, sizePixelsDegs, sigmaPixelsDegs] = smoothingKernelFromConeAperture(...
    theConeApertureDiameterDegs, theRFmapSpatialSupportDegs)

    % Smoothing kernel params
    dxDegs = theRFmapSpatialSupportDegs(2)-theRFmapSpatialSupportDegs(1);
    sigmaDegs = 0.1*theConeApertureDiameterDegs;
    sigmaPixels = sigmaDegs/dxDegs;
    sizePixels = round(6*sigmaDegs/dxDegs);
    if (mod(sizePixels,2) == 0)
        % make it odd
        sizePixels = sizePixels + 1;
    end
    sizePixelsDegs = sizePixels * dxDegs;
    sigmaPixelsDegs = sigmaPixels * dxDegs;
end



function [cLUToriginal, cLUTreversed] = generateReidShapleyRFmapLUT()
    % Reid & Shapley (2002) colormap
    cLUTsampled = [
        155 164 183; ...
        148 156 178; ...
        135 143 170; ...
        121 131 164; ...
        105 115 154; ...
        84 97 142; ...
        52 70 123; ...
        46 57 113 ; ...
        0 0 0; ...
        105 45 45; ...
        132 47 47; ...
        150 49 47; ...
        166 58 58; ...
        175 85 83; ...
        184 112 109; ...
        190 129 129; ...
        195 147 144];

    cLUTsampled = cLUTsampled/max(cLUTsampled(:));

    cLUToriginal(:,1) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,1), linspace(1, size(cLUTsampled,1), 1025));
    cLUToriginal(:,2) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,2), linspace(1, size(cLUTsampled,1), 1025));
    cLUToriginal(:,3) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,3), linspace(1, size(cLUTsampled,1), 1025));
   

    cLUTsampled = [
        46 57 105; ...
        52 70 132; ...
        84 97 150; ...
        105 45 166; ...
        121 131 175; ...
        135 143 184; ...
        148 156 190; ...
        155 164 195; ...
        175 175 175; ...
        195 164 155; ...
        190 156 148; ...
        184 143 135; ...
        175 131 121; ...
        166 45 105; ...
        150 97 84; ...
        132 70 52; ...
        105 57 46];

    cLUTsampled = cLUTsampled/max(cLUTsampled(:));

    cLUTreversed(:,1) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,1), linspace(1, size(cLUTsampled,1), 1025));
    cLUTreversed(:,2) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,2), linspace(1, size(cLUTsampled,1), 1025));
    cLUTreversed(:,3) = interp1(1:size(cLUTsampled,1), cLUTsampled(:,3), linspace(1, size(cLUTsampled,1), 1025));
end


function theRoRincRatio = RoRincrementRatio(theRFmap)
    Ro = sum(theRFmap(:));
    idxNegative = find(theRFmap<0);
    theAllPositiveRFmap = theRFmap;
    theAllPositiveRFmap(idxNegative) = 0;
    Rinc = sum(theAllPositiveRFmap(:));
    theRoRincRatio = Ro/Rinc;
end

function [theRFmap, theRFmapCenter, theRFmapOrientation] = ...
    verticallyAlignedElongationAxisRFmap(theRFmap, theRFmapCenter, theRFmapOrientation)

    if (isempty(theRFmapCenter)) || (isempty(theRFmapOrientation))
        idxNegative = find(theRFmap<0);
        theAllPositiveRFmap = theRFmap;
        theAllPositiveRFmap(idxNegative) = 0;
        theAllPositiveRFmap = imbinarize(theAllPositiveRFmap);
        stats = regionprops(theAllPositiveRFmap, 'image', 'Orientation', 'Centroid');
        if (isempty(theRFmapCenter))
            theRFmapCenter = stats.Centroid - 0.5*size(theRFmap);
        end
        if (isempty(theRFmapOrientation))
            theRFmapOrientation = stats.Orientation+90;
        end
    end

    theRFmap = imtranslate(theRFmap, -theRFmapCenter);
    theRFmap = imrotate(theRFmap, -theRFmapOrientation, 'bilinear','crop');
    theRFmap = imtranslate(theRFmap, theRFmapCenter);
end


function [theAchromaticSlice, theLconeIsolatingSlice, theMconeIsolatingSlice, ...
         theAchromaticPeakRow, theAchromaticPeakCol, ...
         theLconeIsolatingPeakRow, theLconeIsolatingPeakCol] = ...
                    radialyAveragedSlices(theAchromaticRFmap, theLconeIsolatingRFmap, theMconeIsolatingRFmap, ...
                        theRFcenterConeDominance)
    
    % Circularly average the achromatic RF map
    [theAchromaticSlice, theAchromaticPeakCol, theAchromaticPeakRow] = circularlySymmetricSlice(theAchromaticRFmap, []);

    % Circular averaging center
    theForcedPeakRowPeakCol = [theAchromaticPeakCol theAchromaticPeakRow];

    [theLconeIsolatingSlice, theLconeIsolatingPeakCol, theLconeIsolatingPeakRow] = ...
        circularlySymmetricSlice(theLconeIsolatingRFmap, theForcedPeakRowPeakCol);

    [theMconeIsolatingSlice, theMconeIsolatingPeakCol, theMconeIsolatingPeakRow] = ...
        circularlySymmetricSlice(theMconeIsolatingRFmap, theForcedPeakRowPeakCol);

end

function [theSlice, peakCol, peakRow] = circularlySymmetricSlice(theMap, theForcedPeakRowPeakCol)
    [theCircuarlySymmetricMap, peakCol, peakRow] = circularlyAverageRFmap(theMap, theForcedPeakRowPeakCol);
    % Find the center of theCircuarlySymmetricMap
    %[~,idx] = max(theCircuarlySymmetricMap(:));
    %[peakRow, peakCol] = ind2sub(size(theCircuarlySymmetricMap), idx);
    theSlice = theCircuarlySymmetricMap(peakRow,:);
end


function [outRFmap, peakCol, peakRow] = circularlyAverageRFmap(inRFmap, forcedCenter)
    % Circularly average the input RFmap
    % Define quantization. Four was used in early code, but 1 makes more sense.
    quantizationFactor = 1;

    % Make a circularly symmetric version of average optics.
    [m, n] = size(inRFmap);
    if (n ~= m)
        error('RFmap must be a square matrix'); 
    end
    nLinearPixels = m;

    if (isempty(forcedCenter))
        [peakRow, peakCol] = psfFindPeak(inRFmap);
    else
        peakCol = forcedCenter(1);
        peakRow = forcedCenter(2);
    end

    radiusMat = MakeRadiusMat(nLinearPixels, nLinearPixels, peakCol, peakRow);
    outRFmap = zeros(nLinearPixels, nLinearPixels);
    nBands = round(nLinearPixels / quantizationFactor);
    radii = linspace(0, 0.75 * nLinearPixels, nBands);
    for q = 1:length(radii) - 1
        index = find(radiusMat >= radii(q) & radiusMat < radii(q + 1));
        if (~isempty(index))
            outRFmap(index) = mean(inRFmap(index)); 
        end
    end
    outRFmap = sum(inRFmap(:)) * outRFmap / sum(outRFmap(:));
end

