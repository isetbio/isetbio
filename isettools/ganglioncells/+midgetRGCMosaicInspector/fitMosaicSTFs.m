function fitMosaicSTFs(mosaicCenterParams, rfModelParams, opticsParams,  maxRGCsNum)

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);
    
    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');


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
    theMeridianRadius = 1.5 * sqrt(2.0);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', numel(theMeridianAngles), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);

    

    hFigRcDegs = figure(1); clf;
    set(hFigRcDegs, 'Position', [10 10 3000 1250], 'Color', [1 1 1]);
    drawnow;

    hFigRsRcRatios = figure(2); clf;
    set(hFigRsRcRatios, 'Position', [100 100 3000 1250], 'Color', [1 1 1]);
    drawnow;

    hFigSCintSensRatios = figure(3); clf;
    set(hFigSCintSensRatios, 'Position', [200 200 3000 1250], 'Color', [1 1 1]);
    drawnow;

    RGCindices = cell(1, numel(theMeridianAngles));
    signedDistances = cell(1, numel(theMeridianAngles));

    eccentricityLims = 3*[-1 1];
    eccentricityTicks = -5:1:5;

    for iMeridianAngle = 1:numel(theMeridianAngles)
        theMeridianAngle = theMeridianAngles(iMeridianAngle);
        idx = midgetRGCMosaicInspector.rgcIndicesAlongMeridianWithAngle(...
            theMeridianAngle, theMeridianRadius, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs, maxRGCsNum);
        RGCindices{iMeridianAngle} = idx;

        radialDistances = sqrt(sum((theMidgetRGCmosaic.rgcRFpositionsDegs(idx,:)).^2,2));
        if (theMeridianAngle == 90) || (theMeridianAngle == 270)
            signedDistances{iMeridianAngle} = radialDistances  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,2));
        else
            signedDistances{iMeridianAngle} = radialDistances  .* sign(theMidgetRGCmosaic.rgcRFpositionsDegs(idx,1));
        end

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
        ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);

        plotRGCpositionsAlongMeridian(ax1, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesAlongThisMeridian,:), ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));

        RcDegsLims = [0.5 1]; RcDegsTicks = 0.5:0.1:1.0;
        plotDataAlongMeridian(ax2, ax3, ...
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
        ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);

        plotRGCpositionsAlongMeridian(ax1, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesAlongThisMeridian,:), ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));

        RsRcLims = [0 16]; RsRcTicks = 0:2:40;
        plotDataAlongMeridian(ax2, ax3, ...
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
        ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);


        plotRGCpositionsAlongMeridian(ax1, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesAlongThisMeridian,:), ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', numel(rgcIndicesAlongThisMeridian), theMeridianAngle));

        intSensitivitySCLims = [0 1]; intSensitivitySCTicks = 0:0.1:1.0;
        plotDataAlongMeridian(ax2, ax3, ...
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

    hFig =  figure(hFigRcDegs);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_RcDegs.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    hFig =  figure(hFigRsRcRatios);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_RcRsRatios.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    hFig =  figure(hFigSCintSensRatios);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_SCintSensRatios.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    % Append the fittedSTFs structs.
    % Here we use  matfile which lets you read and write to part of variables in a mat file.
    % In that way, if 'theMeridianAngles', 'theMeridianFits' do not exist
    % it will write them to the file, whereas if they do exist it will
    % replace them with the new values
    m = matfile(responsesFileName, 'Writable', true);
    m.theMeridianAngles = theMeridianAngles;
    m.theMeridianFits = theMeridianFits;
    fprintf('Appended theMeridianFits to the responses file: ''%s''.', responsesFileName);
    clear 'm';
end


function plotRGCpositionsAlongMeridian(ax, rgcRFpositionsDegs, ...
    eccentricityLims, eccentricityTicks, plotTitle)

    plot(ax,rgcRFpositionsDegs(idx,1),  rgcRFpositionsDegs(idx,2), 'k.', ...
         'MarkerSize', 15)
    axis(ax1, 'square');
    
    grid(ax1, 'on');
    set(ax1, 'XLim', eccentricityLims, 'YLim', eccentricityLims, ...
            'XTick', eccentricityTicks, 'YTick', eccentricityTicks, ...
            'FontSize', 16);
    title(ax, plotTitle);
    drawnow;
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

    % Allocate memory
    temporalEquivalentEccDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRcDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRsToRcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedKsToKcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedSCintSensRatios = zeros(1, numel(rgcIndicesToAnalyze)); 
    conesNumPooledByTheRFcenter = zeros(1, numel(rgcIndicesToAnalyze)); 
    majorityCenterConeType = zeros(1, numel(rgcIndicesToAnalyze));
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    parfor iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        fprintf('Fitting RGC %d of %d\n', iRGC, numel(rgcIndicesToAnalyze));

        % Do the fitting
        [fittedSTFs{iRGC}, ...
         temporalEquivalentEccDegs(iRGC), ...
         conesNumPooledByTheRFcenter(iRGC), ...
         majorityCenterConeType(iRGC), ...
         achievedRcDegs(iRGC), ...
         achievedRsToRcRatios(iRGC), ...
         achievedKsToKcRatios(iRGC), ...
         achievedSCintSensRatios(iRGC)] = fitSingleRGC(...
                theMidgetRGCmosaic, ...
                theMidgetRGCMosaicResponses, ...
                theRGCindex, ...
                spatialFrequenciesTested, ...
                orientationsTested);

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


function [fittedSTF, temporalEquivalentEccDegs, conesNumPooledByTheRFcenter, majorityCenterConeType, ...
          achievedRcDegs, achievedRsToRcRatio, achievedKsToKcRatio, achievedSCintSensRatio] = ...
            fitSingleRGC(theMidgetRGCmosaic, theMidgetRGCMosaicResponses, theRGCindex, ...
            spatialFrequenciesTested, orientationsTested)
   
        temporalEquivEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:));
        temporalEquivalentEccDegs = sqrt(sum(temporalEquivEccDegs.^2,2));

        % Obtain the responses
        theSingleMidgetRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theRGCindex));

        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = highestExtensionSTF(theSingleMidgetRGCResponses, spatialFrequenciesTested, orientationsTested);

        % Fit the DoG model to the measured STF
        multiStartsNum = 128;
        [RcDegsEstimate, conesNumPooledByTheRFcenter, majorityCenterConeType] = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex);

        [theFittedDoGmodelParams, theDoGmodelFitOfTheMeasuredSTF] = ...
            RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTF, ...
                    RcDegsEstimate, ...
                    multiStartsNum);


        % The RcDegs
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RcDegs'));
        achievedRcDegs = theFittedDoGmodelParams.finalValues(idx);

        % The Rs/Rc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RsToRc'));
        achievedRsToRcRatio = theFittedDoGmodelParams.finalValues(idx);

        % The Ks/Kc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'kS/kC'));
        achievedKsToKcRatio = theFittedDoGmodelParams.finalValues(idx);

        % The S/C int sens ratio
        achievedSCintSensRatio = achievedKsToKcRatio * (achievedRsToRcRatio)^2;


        % Determine whether the majority of center cones are L- or M-
        theMajorityCenterConeType = theMidgetRGCmosaic.majorityCenterConeType(theRGCindex);
        
        % Compute the indices of the triangulating RTVFobjects and their
        % contributing weights

        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theRGCindex);

        pipelineScaleFactorBasedOnLowestSF = 0;
        triangulatingRTVFobjSTFdata = cell(1, numel(triangulatingRTVFobjIndices));

        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
            iObj  = triangulatingRTVFobjIndices(iNearbyObj);
            theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{iObj};

            switch (theMajorityCenterConeType)
                case cMosaic.LCONE_ID
                    theRFcomputeStruct = theRTVFTobj.LconeRFcomputeStruct;

                case cMosaic.MCONE_ID
                    theRFcomputeStruct = theRTVFTobj.MconeRFcomputeStruct;

                otherwise
                    error('How can the majority cone type be not L- or M- ??')
            end

            % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
            scaleFactorBasedOnLowestSF = theRFcomputeStruct.theSTF.target(1)/theRFcomputeStruct.theSTF.fitted(1);
            pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF + theRFcomputeStruct.theSTF.target(1)*triangulatingRTVFobjWeights(iNearbyObj);
            
            % The spatial support
            triangulatingRTVFobjSTFdata{iNearbyObj} = struct();
            triangulatingRTVFobjSTFdata{iNearbyObj}.spatialFrequencySupport = theRFcomputeStruct.theSTF.support(:);
            % The model-achieved STF
            triangulatingRTVFobjSTFdata{iNearbyObj}.fittedSTF = theRFcomputeStruct.theSTF.fitted(:)*scaleFactorBasedOnLowestSF;
            % The model-target STF
            triangulatingRTVFobjSTFdata{iNearbyObj}.targetSTF = theRFcomputeStruct.theSTF.target(:);
        end

        
        fittedSTF = struct(...
            'targetRGC', theRGCindex, ...
            'theMultiFocalRTVFmodelSTFdata', triangulatingRTVFobjSTFdata, ...
            'theMultiFocalRTVFmodelWeights', triangulatingRTVFobjWeights, ...
            'spatialFrequencySupport', spatialFrequenciesTested, ...
            'theMeasuredSTF', theMeasuredSTF, ... 
            'allMeasuredSTFs', allMeasuredSTFs, ...
            'theDoGmodelFitOfTheMeasuredSTF', theDoGmodelFitOfTheMeasuredSTF, ...
            'theFittedDoGmodelParams', theFittedDoGmodelParams);
            
    
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

