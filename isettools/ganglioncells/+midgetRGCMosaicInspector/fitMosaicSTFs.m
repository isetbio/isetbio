function fitMosaicSTFs(mosaicCenterParams, mosaicSurroundParams,  maxRGCsNum)

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, mosaicSurroundParams.H1cellIndex);
    
    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');

    % RGC indices to fit
%     rgcIndicesToAnalyze = midgetRGCMosaicInspector.rgcIndicesToAnalyze(...
%         frozenMosaicFileName, ...
%         'maxRGCsNum', maxRGCsNum);

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


    % Examine every 15 degs
    deltaAngle = 22.5;
    theMeridianAngles  = 0:deltaAngle:(180-deltaAngle);
    radius = 1.5 * sqrt(2.0);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', numel(theMeridianAngles), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);

    

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 3000 850], 'Color', [1 1 1]);
    drawnow;

    hFigRcDegs = figure(2); clf;
    set(hFigRcDegs, 'Position', [10 10 3000 850], 'Color', [1 1 1]);
    drawnow;

    hFigRsRcRatios = figure(3); clf;
    set(hFigRsRcRatios, 'Position', [100 100 3000 850], 'Color', [1 1 1]);
    drawnow;

    hFigSCintSensRatios = figure(4); clf;
    set(hFigSCintSensRatios, 'Position', [200 200 3000 850], 'Color', [1 1 1]);
    drawnow;

    RGCindices = cell(1, numel(theMeridianAngles));
    signedDistances = cell(1, numel(theMeridianAngles));

    for iMeridianAngle = 1:numel(theMeridianAngles)
        theMeridianAngle = theMeridianAngles (iMeridianAngle);
        x = radius * cosd(theMeridianAngle );
        y = radius * sind(theMeridianAngle );
    
        theROI = regionOfInterest('shape', 'line', 'from', [x,y], 'to', [-x,-y], 'thickness', 0.01);
        samplingPoints = 400;  % sample the perimeter of the ROI along 1000 points');
        pointsPerSample = 10;  % return up to 10 points for each sample along the perimeter');
        maxDistance = 0.1;     % points must be no further than 0.1 degs away from the closest perimeter sample');
        idx = theROI.indicesOfPointsAround(theMidgetRGCmosaic.rgcRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
        
        if (~isempty(maxRGCsNum))
            skip = floor(numel(idx)/maxRGCsNum);
        else
            skip = 1;
        end
        idx = idx(1:skip:numel(idx));
        RGCindices{iMeridianAngle} = idx;

        radialDistances = sqrt(sum((theMidgetRGCmosaic.rgcRFpositionsDegs(idx,:)).^2,2));
        if (theMeridianAngle == 90)
            signedDistances{iMeridianAngle} = radialDistances  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,2));
        else
            signedDistances{iMeridianAngle} = radialDistances  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,1));
        end

        ax = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
        plot(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,1), ...
             theMidgetRGCmosaic.rgcRFpositionsDegs(idx,2), 'k.', ...
             'MarkerSize', 15)
        axis(ax, 'square');
        eccentricityLims = 3*[-1 1];
        eccentricityTicks = -5:1:5;
        grid(ax, 'on');
        set(ax, 'XLim', eccentricityLims, 'YLim', eccentricityLims, ...
                'XTick', eccentricityTicks, 'YTick', eccentricityTicks, ...
                'FontSize', 16);
        title(ax, sprintf('%d RGCs along the\n%2.1f deg meridian', numel(idx), theMeridianAngle));
        drawnow;
    end

    

    % Retrieve the C&K data
    [CronerKaplanRcDegs.temporalEccDegs, CronerKaplanRcDegs.val] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    [CronerKaplanRcRsRatios.temporalEccDegs, CronerKaplanRcRsRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

    [CronerKaplanSCintSensRatios.temporalEccDegs, CronerKaplanSCintSensRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();

    [CronerKaplanKsKcRatios.temporalEccDegs, CronerKaplanKsKcRatios.val] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

    theMeridianFits = cell(1, numel(theMeridianAngles));
    for iMeridianAngle = 1:numel(theMeridianAngles)

        theMeridianAngle = theMeridianAngles(iMeridianAngle);
        
        % Fit all RGCs located at this meridian
        rgcIndicesAlongThisMeridian = RGCindices{iMeridianAngle};
        [fittedParams, fittedSTFs] = fitSelectSTFs(rgcIndicesAlongThisMeridian, ...
            spatialFrequenciesTested, orientationsTested, ...
            theMidgetRGCmosaic, theMidgetRGCMosaicResponses);

        theMeridianFits{iMeridianAngle} = struct(...
            'rgcIndicesAlongThisMeridian', rgcIndicesAlongThisMeridian, ...
            'fittedParams', fittedParams, ...
            'fittedSTFs', fittedSTFs);



        % The Rc degs values
        figure(hFigRcDegs);
        ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
        ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);

        RcDegsLims = [0.5 1]; RcDegsTicks = 0.5:0.1:1.0;
        plotDataAlongMeridian(ax1, ax2, ...
            signedDistances{iMeridianAngle}, fittedParams.achievedRcDegs*60, ...
            eccentricityLims, eccentricityTicks, ...
            RcDegsLims, RcDegsTicks, ...
            [RcDegsLims(1) 8], 0.5:0.5:10, ...
            'RcDegs (arc min)', sprintf('%d RGCs along %2.1f degs', numel(rgcIndicesAlongThisMeridian), theMeridianAngle), ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.majorityCenterConeType, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanRcDegs.temporalEccDegs, CronerKaplanRcDegs.val*60 ...
            );
        

        % The Rs/Rc ratio values
        figure(hFigRsRcRatios);
        ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
        ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);

        RsRcLims = [0 16]; RsRcTicks = 0:2:40;
        plotDataAlongMeridian(ax1, ax2, ...
            signedDistances{iMeridianAngle}, fittedParams.achievedRsToRcRatios, ...
            eccentricityLims, eccentricityTicks, ...
            RsRcLims, RsRcTicks, ...
            [RsRcLims(1) 32], RsRcTicks, ...
            'Rs/Rc ratio', sprintf('%d RGCs along %2.1f degs', numel(rgcIndicesAlongThisMeridian), theMeridianAngle), ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.majorityCenterConeType, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanRcRsRatios.temporalEccDegs, 1./CronerKaplanRcRsRatios.val...
            );

        % The S/C integrated sensitivity ratios
        figure(hFigSCintSensRatios);
        ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
        ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);

        intSensitivitySCLims = [0 1]; intSensitivitySCTicks = 0:0.1:1.0;
        plotDataAlongMeridian(ax1, ax2, ...
            signedDistances{iMeridianAngle}, fittedParams.achievedSCintSensRatios, ...
            eccentricityLims, eccentricityTicks, ...
            intSensitivitySCLims , intSensitivitySCTicks , ...
            intSensitivitySCLims , intSensitivitySCTicks , ...
            'S/C int. sens. ratio', sprintf('%d RGCs along\nthe %2.1f degs meridian', numel(rgcIndicesAlongThisMeridian), theMeridianAngle), ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.majorityCenterConeType, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanSCintSensRatios.temporalEccDegs, CronerKaplanSCintSensRatios.val ...
            );

    end

    % Append the fittedSTFs structs
    save(responsesFileName, 'theMeridianFits', 'theMeridianAngles', '-append');
end

function plotDataAlongMeridian(ax1, ax2, ...
            signedDistances, fittedParams, ...
            XLims, XTicks, ...
            YLims, YTicks, ...
            YLimsCK, YTicksCK, ...
            yLabelString, titleString, ...
            conesNumPooledByTheRFcenters, ...
            majorityCenterConeTypes, ...
            temporalEquivalentEccDegs, ...
            CronerKaplanTemporalEccDegs, CronerKaplanValues)
 

        % Plot on ax1
        majorityCenterConeType = unique(majorityCenterConeTypes);
        for i = 1:numel(majorityCenterConeType)
             rgcIndicesToPlot = find(majorityCenterConeTypes == majorityCenterConeType(i));
             switch (majorityCenterConeType(i))
                case cMosaic.LCONE_ID
                    markerColor = [1 0 0];
                case cMosaic.MCONE_ID
                    markerColor = [0 0.7 0];
                otherwise
                    markerColor = [0 0 0];
            end
            plot(ax1, signedDistances(rgcIndicesToPlot), fittedParams(rgcIndicesToPlot), '.', 'MarkerSize', 15, ...
                'MarkerFaceColor', markerColor, 'MarkerEdgeColor', markerColor);
            hold(ax1, 'on');
        end
        
        axis(ax1, 'square');
        set(ax1, 'XLim', XLims, 'XTick', XTicks, 'YLim', YLims, 'YTick', YTicks, 'FontSize', 16);
        grid(ax1, 'on');
        xtickangle(ax1, 0);
        ylabel(ax1, yLabelString);
        xlabel(ax1, 'eccentricity (degs)');
        title(ax1, titleString);

        % Now replot on ax2, along with the Croner & Kaplan data
        for i = 1:numel(majorityCenterConeType)
             rgcIndicesToPlot = find(majorityCenterConeTypes == majorityCenterConeType(i));
             switch (majorityCenterConeType(i))
                case cMosaic.LCONE_ID
                    markerColor = [1 0 0];
                case cMosaic.MCONE_ID
                    markerColor = [0 0.7 0];
                otherwise
                    markerColor = [0 0 0];
            end
            plot(ax2, temporalEquivalentEccDegs(rgcIndicesToPlot), fittedParams(rgcIndicesToPlot), '.', 'MarkerSize', 15, ...
                'MarkerFaceColor', markerColor, 'MarkerEdgeColor', markerColor);
            hold(ax2, 'on');
        end

        scatter(ax2, CronerKaplanTemporalEccDegs, CronerKaplanValues, 121, 'o', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 1.0);
        axis(ax2, 'square');
        set(ax2, 'XLim', [0.01 30], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], 'XScale', 'log');
        xtickangle(ax2, 0);
        set(ax2, 'YLim', YLimsCK, 'YTick', YTicksCK, 'FontSize', 16);
        grid(ax2, 'on');
        ylabel(ax2, yLabelString);
        xlabel(ax2, 'temporal equiv. eccentricity (degs)');
        drawnow

end


function [d, fittedSTFs] = fitSelectSTFs(rgcIndicesToAnalyze, ...
    spatialFrequenciesTested, orientationsTested, ...
    theMidgetRGCmosaic, theMidgetRGCMosaicResponses)

    temporalEquivalentEccDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRcDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRsToRcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedKsToKcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedSCintSensRatios = zeros(1, numel(rgcIndicesToAnalyze)); 
    conesNumPooledByTheRFcenter = zeros(1, numel(rgcIndicesToAnalyze)); 
    majorityCenterConeType = zeros(1, numel(rgcIndicesToAnalyze)); 

    % Allocate memory
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        fprintf('Fitting RGC %d of %d, located at (%2.2f,%2.2f degs)\n', ...
            iRGC, numel(rgcIndicesToAnalyze), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,1), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,2));

        temporalEquivEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:));
        temporalEquivalentEccDegs(iRGC) = sqrt(sum(temporalEquivEccDegs.^2,2));

        % Obtain the responses
        theSingleMidgetRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theRGCindex));

        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = highestExtensionSTF(theSingleMidgetRGCResponses, spatialFrequenciesTested, orientationsTested);

        % Fit the DoG model to the measured STF
        multiStartsNum = 128;
        [RcDegsEstimate, conesNumPooledByTheRFcenter(iRGC), majorityCenterConeType(iRGC)] = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex);

        [theFittedDoGmodelParams, theDoGmodelFitOfTheMeasuredSTF] = ...
            RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTF, ...
                    RcDegsEstimate, ...
                    multiStartsNum);


        % The RcDegs
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RcDegs'));
        achievedRcDegs(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The Rs/Rc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RsToRc'));
        achievedRsToRcRatios(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The Ks/Kc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'kS/kC'));
        achievedKsToKcRatios(iRGC) = theFittedDoGmodelParams.finalValues(idx);

        % The S/C int sens ratio
        achievedSCintSensRatios(iRGC) = achievedKsToKcRatios(iRGC) * (achievedRsToRcRatios(iRGC))^2;


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

          

            subplot(2,4,[5 6]);
            scatter(temporalEquivalentEccDegs, achievedRsToRcRatios, 140, ...
                'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeAlpha', 0.5);
    
            set(gca, 'XScale', 'log', 'XLim', [0.01 3], 'XTick', [0.01 0.03 0.1 0.3 1 3]);
            set(gca, 'YLim', [0 16], 'YTick', 0:2:16);
            set(gca, 'FontSize', 16);
            grid on
            xlabel('temporal eccentricity (degs)');
            ylabel('Rs/Rc ratio')
    
    
            subplot(2,4,[7 8]);
            scatter(temporalEquivalentEccDegs, achievedSCintSensRatios, 140, ...
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

    
    d = struct(...
        'temporalEquivalentEccDegs', temporalEquivalentEccDegs, ...
        'conesNumPooledByTheRFcenter', conesNumPooledByTheRFcenter, ...
        'majorityCenterConeType', majorityCenterConeType, ...
        'achievedRcDegs', achievedRcDegs, ...
        'achievedRsToRcRatios', achievedRsToRcRatios, ...
        'achievedKsToKcRatios', achievedKsToKcRatios, ...
        'achievedSCintSensRatios', achievedSCintSensRatios ...
        );

end

function [RcDegs, conesNumPooledByTheRFcenter,  majorityCenterConeType] = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex)
    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);
    centerConeTypes = theMidgetRGCmosaic.inputConeMosaic.coneTypes(indicesOfCenterCones);
    lConesNum = numel(find(centerConeTypes == cMosaic.LCONE_ID));
    mConesNum = numel(find(centerConeTypes == cMosaic.MCONE_ID));
    if (lConesNum > mConesNum)
        majorityCenterConeType = cMosaic.LCONE_ID;
    elseif (lConesNum < mConesNum)
        majorityCenterConeType = cMosaic.MCONE_ID;
    else
        majorityCenterConeType = 0;
    end

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

