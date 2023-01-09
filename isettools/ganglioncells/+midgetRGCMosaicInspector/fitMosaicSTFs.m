function fitMosaicSTFs(mosaicCenterParams, mosaicSurroundParams,  maxRGCsNum)

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, mosaicSurroundParams.H1cellIndex);
    
    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');

    % RGC indices to fit
    rgcIndicesToAnalyze = midgetRGCMosaicInspector.rgcIndicesToAnalyze(...
        frozenMosaicFileName, ...
        'maxRGCsNum', maxRGCsNum);

    % Ask the user to specify the optics position for which responses were saved
    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        opticsPositionDegs = input('\nEnter the optics position that was used to compute the responses ([x y]): ');
    end

    % Generate the responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);
       

    % Load the mosaic responses
    load(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',  0.08, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.05);

    % Examine every 22.5 degs
    theMeridianAngles  = 0:22.5:(180-22.5);
    radius = 1.5 * sqrt(2.0);

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1400 700]);

    for iMeridianAngle = 1:numel(theMeridianAngles )
        theMeridianAngle = theMeridianAngles (iMeridianAngle);
        x = radius * cosd(theMeridianAngle );
        y = radius * sind(theMeridianAngle );
    
        theROI = regionOfInterest('shape', 'line', 'from', [x,y], 'to', [-x,-y], 'thickness', 0.01);
        samplingPoints = 200;  % sample the perimeter of the ROI along 1000 points');
        pointsPerSample = 10;  % return up to 10 points for each sample along the perimeter');
        maxDistance = 0.1;     % points must be no further than 0.1 degs away from the closest perimeter sample');
        idx = theROI.indicesOfPointsAround(theMidgetRGCmosaic.rgcRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
        
        radialDistance = sqrt(sum((theMidgetRGCmosaic.rgcRFpositionsDegs(idx,:)).^2,2));

        if (theMeridianAngle == 90)
            signedDistances{iMeridianAngle} = radialDistance  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,2));
        else
            signedDistances{iMeridianAngle} = radialDistance  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,1));
        end
        RGCindices{iMeridianAngle} = idx;

        i = floor((iMeridianAngle-1)/4)+1;
        j = mod(iMeridianAngle-1,4)+1;
        ax = subplot('Position', subplotPosVectors(i,j).v);
        plot(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,1), ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(idx,2), 'k.');
        axis(ax, 'square');
        set(ax, 'XLim', 1.5*[-1 1], 'YLim', 1.5*[-1 1]);
        title(ax, sprintf('%d RGCs along %2.1f degs', numel(idx), theMeridianAngle ));
    end

    hFigRcDegs = figure(2); clf;
    set(hFigRcDegs, 'Position', [10 10 1400 700]);
    drawnow;

    hFigRsRcRatios = figure(3); clf;
    set(hFigRsRcRatios, 'Position', [100 100 1400 700]);
    drawnow;

    hFigSCintSensRatios = figure(4); clf;
    set(hFigSCintSensRatios, 'Position', [200 200 1400 700]);
    drawnow;

    theMeridianFits = cell(1, numel(theMeridianAngles));
    for iMeridianAngle = 1:numel(theMeridianAngles)

        theMeridianAngle = theMeridianAngles(iMeridianAngle);
        i = floor((iMeridianAngle-1)/4)+1;
        j = mod(iMeridianAngle-1,4)+1;
        
        % Fit all RGCs located at this meridian
        rgcIndicesAlongThisMeridian = RGCindices{iMeridianAngle};
        [fittedParams, fittedSTFs] = fitSelectSTFs(rgcIndicesAlongThisMeridian, ...
            spatialFrequenciesTested, orientationsTested, ...
            theMidgetRGCmosaic, theMidgetRGCMosaicResponses);

        theMeridianFits{iMeridianAngle} = struct(...
            'rgcIndicesAlongThisMeridian', rgcIndicesAlongThisMeridian, ...
            'fittedParams', fittedParams, ...
            'fittedSTFs', fittedSTFs);


        figure(hFigRcDegs);
        ax = subplot('Position', subplotPosVectors(i,j).v);
        plot(ax, signedDistances{iMeridianAngle}, fittedParams.achievedRcDegs, 'r.');
        axis(ax, 'square');
        set(ax, 'XLim', 1.5*sqrt(2)*[-1 1]);
        ylabel(ax, 'RcDegs');
        xlabel(ax, 'eccentricity (degs)');
        title(ax, sprintf('%d RGCs along %2.1f degs', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));
        drawnow

        figure(hFigRsRcRatios);
        ax = subplot('Position', subplotPosVectors(i,j).v);
        plot(ax, signedDistances{iMeridianAngle}, fittedParams.achievedRsToRcRatios, 'r.');
        axis(ax, 'square');
        set(ax, 'XLim', 1.5*sqrt(2)*[-1 1], 'YLim', [0 16], 'YTick', 0:2:16);
        ylabel(ax, 'Rs/Rc ratio');
        xlabel(ax, 'eccentricity (degs)');
        title(ax, sprintf('%d RGCs along %2.1f degs', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));
        drawnow

        figure(hFigSCintSensRatios);
        ax = subplot('Position', subplotPosVectors(i,j).v);
        plot(ax, signedDistances{iMeridianAngle}, fittedParams.achievedSCintSensRatios, 'r.');
        axis(ax, 'square');
        set(ax, 'XLim', 1.5*sqrt(2)*[-1 1], 'YLim', [0 1], 'YTick', 0:0.2:1.0);
        ylabel(ax, 'S/C int. sens. ratio');
        xlabel(ax, 'eccentricity (degs)');
        title(ax, sprintf('%d RGCs along %2.1f degs', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));
        drawnow
    end

    % Append the fittedSTFs structs
    save(responsesFileName, 'theMeridianFits', 'theMeridianAngles', '-append');


end

function [d, fittedSTFs] = fitSelectSTFs(rgcIndicesToAnalyze, ...
    spatialFrequenciesTested, orientationsTested, ...
    theMidgetRGCmosaic, theMidgetRGCMosaicResponses)

    d = struct(...
        'temporalEquivalentEccDegs', zeros(1, numel(rgcIndicesToAnalyze)), ...
        'achievedRcDegs', zeros(1, numel(rgcIndicesToAnalyze)), ...
        'achievedRsToRcRatios', zeros(1, numel(rgcIndicesToAnalyze)), ...
        'achievedKsToKcRatios', zeros(1, numel(rgcIndicesToAnalyze)), ...
        'achievedSCintSensRatios', zeros(1, numel(rgcIndicesToAnalyze)) ...
        );

    % Allocate memory
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        fprintf('Fitting RGC %d of %d, located at (%2.2f,%2.2f degs)\n', ...
            iRGC, numel(rgcIndicesToAnalyze), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,1), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,2));

        temporalEquivEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:));
        d.temporalEquivalentEccDegs(iRGC) = sqrt(sum(temporalEquivEccDegs.^2,2));

        % Obtain the responses
        theSingleMidgetRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theRGCindex));

        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = highestExtensionSTF(theSingleMidgetRGCResponses, spatialFrequenciesTested, orientationsTested);

        % Fit the DoG model to the measured STF
        multiStartsNum = 128;
        [theFittedDoGmodelParams, theDoGmodelFitOfTheMeasuredSTF] = ...
            RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTF, ...
                    retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex), ...
                    multiStartsNum);


        % The RcDegs
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RcDegs'));
        d.achievedRcDegs(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The Rs/Rc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RsToRc'));
        d.achievedRsToRcRatios(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The Ks/Kc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'kS/kC'));
        d.achievedKsToKcRatios(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The S/C int sens ratio
        d.achievedSCintSensRatios(iRGC) = d.achievedKsToKcRatios(iRGC) * (d.achievedRsToRcRatios(iRGC))^2;


        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theRGCindex);

        pipelineScaleFactorBasedOnLowestSF = 0;
        triangulatingRTVFobjSTFdata = cell(1, numel(triangulatingRTVFobjIndices));

        for iObj = 1:numel(triangulatingRTVFobjIndices)
            theRTVFobjIndex = triangulatingRTVFobjIndices(iObj);
            theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{theRTVFobjIndex};

            % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
            scaleFactorBasedOnLowestSF = theRTVFTobj.rfComputeStruct.theSTF.target(1)/theRTVFTobj.rfComputeStruct.theSTF.fitted(1);
            pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF + theRTVFTobj.rfComputeStruct.theSTF.target(1)*triangulatingRTVFobjWeights(iObj);
            
            % The spatial support
            triangulatingRTVFobjSTFdata{iObj} = struct();
            triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport = theRTVFTobj.rfComputeStruct.theSTF.support(:);
            % The model-achieved STF
            triangulatingRTVFobjSTFdata{iObj}.fittedSTF = theRTVFTobj.rfComputeStruct.theSTF.fitted(:)*scaleFactorBasedOnLowestSF;
            % The model-target STF
            triangulatingRTVFobjSTFdata{iObj}.targetSTF = theRTVFTobj.rfComputeStruct.theSTF.target(:);
        end

        
    
        visualizeSTFfits = false;
        if (visualizeSTFfits)

            % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
            pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF /theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes(1);

            hFig = figure(10000); clf;
            set(hFig, 'Position',[10 10 1800 1100]);
    
            for iObj = 1:numel(triangulatingRTVFobjIndices) 
    
                RTVFobjPosition = theMidgetRGCmosaic.theSamplingPositionGrid(triangulatingRTVFobjIndices(iObj),:);
    
                subplot(2,4, iObj)
                plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.targetSTF, ...
                    'k-', 'LineWidth', 1.5);
                hold on;
                plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.fittedSTF, ...
                    'r-', 'LineWidth', 1.5);
                plot(spatialFrequenciesTested, theMeasuredSTF*pipelineScaleFactorBasedOnLowestSF, 'co-', 'MarkerSize', 14, 'LineWidth', 1.5, 'MarkerFaceColor', [0 1 1], 'Color', [0 0.6 0.6]);
                plot(theDoGmodelFitOfTheMeasuredSTF.sfHiRes, theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes*pipelineScaleFactorBasedOnLowestSF, 'b-', 'LineWidth', 1.5);
                axis 'square'
                legend({'RTVF: target', 'RTVF: fitted', 'pipeline: measured', 'pipeline: DoGmodel'}, 'Location', 'SouthWest');
                set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', ...
                    'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.1 100]);
                grid on
                set(gca, 'FontSize', 16);
                xlabel('spatial frequency (c/deg)')
                title(sprintf('component RTVF model (#%d)\nweight %2.3f, located at (%2.2f,%2.2f) degs', ...
                    triangulatingRTVFobjIndices(iObj), triangulatingRTVFobjWeights(iObj), RTVFobjPosition(1), RTVFobjPosition(2)));
            end
        
            subplot(2,4,4);
            plot(spatialFrequenciesTested, theMeasuredSTF*pipelineScaleFactorBasedOnLowestSF, 'co-', 'MarkerSize', 14, 'LineWidth', 1.5, 'MarkerFaceColor', [0 1 1], 'Color', [0 0.6 0.6]);
            hold on
            plot(theDoGmodelFitOfTheMeasuredSTF.sfHiRes, theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes*pipelineScaleFactorBasedOnLowestSF, 'b-', 'LineWidth', 1.5);
            legend({'pipeline: measured', 'pipeline: DoGmodel'}, 'Location', 'SouthWest');
            axis 'square'
            set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.1 100]);
            grid on
            set(gca, 'FontSize', 16);
            xlabel('spatial frequency (c/deg)')
            title(sprintf('RGC %d of %d at (%2.2f,%2.2f degs)', ...
                iRGC, numel(rgcIndicesToAnalyze), ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,1), ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,2)));

            idx = find(strcmp(theFittedDoGmodelParams.names, 'RsToRc'));
            theAchievedRsToRcRatio(iRGC) = theFittedDoGmodelParams.finalValues(idx);
            idx = find(strcmp(theFittedDoGmodelParams.names, 'kS/kC'));
            theAchievedKsToKcRatio = theFittedDoGmodelParams.finalValues(idx);
            theAchievedSCintSensRatio(iRGC) = theAchievedKsToKcRatio *(theAchievedRsToRcRatio(iRGC))^2;

            subplot(2,4,[5 6]);
            scatter(theTemporalEccDegs, theAchievedRsToRcRatio, 140, ...
                'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeAlpha', 0.5);
    
            set(gca, 'XScale', 'log', 'XLim', [0.01 3], 'XTick', [0.01 0.03 0.1 0.3 1 3]);
            set(gca, 'YLim', [0 16], 'YTick', 0:2:16);
            set(gca, 'FontSize', 16);
            grid on
            xlabel('temporal eccentricity (degs)');
            ylabel('Rs/Rc ratio')
    
    
            subplot(2,4,[7 8]);
            scatter(theTemporalEccDegs, theAchievedSCintSensRatio, 140, ...
                'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeAlpha', 0.5);
            set(gca, 'XScale', 'log', 'XLim', [0.01 3], 'XTick', [0.01 0.03 0.1 0.3 1 3]);
            set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1.0);
            xlabel('temporal eccentricity (degs)');
            ylabel('S/C int. sensitivity ratio');
            set(gca, 'FontSize', 16);
            grid on
        end



        fittedSTFs{iRGC} = struct(...
            'targetRGC', theRGCindex, ...
            'theMultiFocalRTVFmodelSTFdata', triangulatingRTVFobjSTFdata, ...
            'theMultiFocalRTVFmodelWeights', triangulatingRTVFobjWeights, ...
            'spatialFrequencySupport', spatialFrequenciesTested, ...
            'theMeasuredSTF', theMeasuredSTF, ... 
            'allMeasuredSTFs', allMeasuredSTFs, ...
            'theDoGmodelFitOfTheMeasuredSTF', theDoGmodelFitOfTheMeasuredSTF, ...
            'theFittedDoGmodelParams', theFittedDoGmodelParams);
            
    end % iRGC
end

function RcDegs = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex)
    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);
    conesNumPooledByTheRFcenter = numel(indicesOfCenterCones);
    coneRcDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones)) * ...
                 theMidgetRGCmosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
    RcDegs = sqrt(conesNumPooledByTheRFcenter)*coneRcDegs;
end

function [theHighestExtensionSTF, theMeasuredSTFs] = highestExtensionSTF(theMidgetRGCMosaicResponses, spatialFrequenciesTested, orientationsTested)

    % Allocate memory
    theMeasuredSTFs = zeros(numel(orientationsTested),numel(spatialFrequenciesTested));

    for iSF = 1:numel(spatialFrequenciesTested)
        for iOri = 1:numel(orientationsTested)
            % Retrieve the mRGC response time-series
            theResponseTimeSeries = squeeze(theMidgetRGCMosaicResponses(iOri, iSF, :));
    
            % Compute the response modulation for this SF
            theMeasuredSTFs(iOri, iSF) = max(theResponseTimeSeries)-min(theResponseTimeSeries);
        end % iORI
    end % iSF

    theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    spatialFrequenciesInterpolated = linspace(spatialFrequenciesTested(1),spatialFrequenciesTested(end), 50);

    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to 20% of max
        theSTFatThisOri = squeeze(theMeasuredSTFs(iOri,:));
        theSTFatThisOriInterpolated = interp1(spatialFrequenciesTested, theSTFatThisOri, spatialFrequenciesInterpolated);
        [mag, iSFpeak] = max(theSTFatThisOri);
        thresholdSTF = mag * 0.2;

        ii = iSFpeak;
        keepGoing = true; iStop = [];
        while (ii < numel(spatialFrequenciesInterpolated)-1)&&(keepGoing)
            ii = ii + 1;
            if (theSTFatThisOriInterpolated(ii)>=thresholdSTF) && (theSTFatThisOriInterpolated(ii+1)<thresholdSTF)
                keepGoing = false;
                iStop = ii;
            end
        end % while
        if (~isempty(iStop))
            maxSF(iOri) = spatialFrequenciesInterpolated(iStop);
        end
    end % iOri

    % Best orientation
    if (any(isnan(maxSF)))
        theSTFatTheHighestSF = squeeze(theMeasuredSTFs(:,end));
        [~, iBestOri] = max(theSTFatTheHighestSF(:));
    else
        [~, iBestOri] = max(maxSF);
    end
    theHighestExtensionSTF = squeeze(theMeasuredSTFs(iBestOri,:));
end

