function fitMosaicSTFs(mosaicCenterParams, mosaicSurroundParams, maxRGCsNum)

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.generateMosaicFileName(...
        mosaicCenterParams);
   
    % Generate responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileNameForMosaicFileName(...
        mosaicFileName, mosaicSurroundParams.H1cellIndex);
 
    % RGC indices to fit
    rgcIndicesToAnalyze = midgetRGCMosaicInspector.rgcIndicesToAnalyze(...
        mosaicFileName, ...
        'maxRGCsNum', maxRGCsNum);

    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the mosaic responses
    load(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');


    % Allocate memory
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        fprintf('Fitting RGC %d of %d, located at (%2.2f,%2.2f degs)\n', ...
            iRGC, numel(rgcIndicesToAnalyze), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,1), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,2));

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

        % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
        pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF /theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes(1);

        figure(10000); clf;
        for iObj = 1:numel(triangulatingRTVFobjIndices) 

            RTVFobjPosition = theMidgetRGCmosaic.theSamplingPositionGrid(triangulatingRTVFobjIndices(iObj),:);

            subplot(1,numel(triangulatingRTVFobjIndices), iObj)
            plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.targetSTF, ...
                'k-', 'LineWidth', 1.5);
            hold on;
            plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.fittedSTF, ...
                'r-', 'LineWidth', 1.5);
            plot(spatialFrequenciesTested, theMeasuredSTF*pipelineScaleFactorBasedOnLowestSF, 'co-', 'MarkerSize', 14, 'LineWidth', 1.5, 'MarkerFaceColor', [0 1 1], 'Color', [0 0.6 0.6]);
            plot(theDoGmodelFitOfTheMeasuredSTF.sfHiRes, theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes*pipelineScaleFactorBasedOnLowestSF, 'b-', 'LineWidth', 1.5);
            legend({'RTVF: target', 'RTVF: fitted', 'pipeline: measured', 'pipeline: DoGmodel'}, 'Location', 'SouthWest');
            set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.01 100]);
            grid on
            set(gca, 'FontSize', 16);
            xlabel('spatial frequency (c/deg)')
            title(sprintf('component RTVF model\nposition:(%2.2f,%2.2f degs), weight: %2.3f', ...
                RTVFobjPosition(1), RTVFobjPosition(2), triangulatingRTVFobjWeights(iObj)));
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

    % Append the fittedSTFs structs
    save(responsesFileName, 'fittedSTFs', 'rgcIndicesToAnalyze', '-append');

end

function RcDegs = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex)
    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
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

