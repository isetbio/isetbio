function compareRTVFmodelSTFtoMeasuredSTF(mosaicCenterParams,  rfModelParams, opticsParams)
    
    % Generate mosaic filename and directory
    mosaicFileName = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');
    theUnfrozenMidgetRGCMosaic = theMidgetRGCmosaic;
    

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);

    load(frozenMosaicFileName, 'theMidgetRGCmosaic');
    theFrozenMidgetRGCMosaic = theMidgetRGCmosaic;
    
    
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(...
        mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    % Load the computed R2VFTobjects
    load(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');



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

    for i = 1:size(theOpticsPositionGrid,1)
        fprintf('RTVF(%d) was fitted at %2.3f,%2.3f\n', ...
            i,theOpticsPositionGrid(i,1), theOpticsPositionGrid(i,2));
    end
    theIndex = input('Enter RTVF obj index to analyze: ');

    [theSourceRGCindex, theCenterConeTypes] = sourceRGCindexForRTVFobj(theRTFVTobjList{theIndex}, ...
            theUnfrozenMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix, ...
            theUnfrozenMidgetRGCMosaic.inputConeMosaic.coneTypes);

    d = sqrt(sum((theUnfrozenMidgetRGCMosaic.rgcRFpositionsDegs).^2,2));
    theSourceRGCindex
    nearbySourceRGCindices = find(d < 0.05)

    analyzePosition(1,...
        theRTFVTobjList{theIndex}, ...
        theSourceRGCindex, ...
        theCenterConeTypes, ...
        theUnfrozenMidgetRGCMosaic, ...
        theMidgetRGCMosaicResponses, ...
        spatialFrequenciesTested, ...
        orientationsTested);

    for i = 1:numel(nearbySourceRGCindices)
        analyzePosition(1+i,...
            theRTFVTobjList{theIndex}, ...
            nearbySourceRGCindices(i), ...
            theCenterConeTypes, ...
            theUnfrozenMidgetRGCMosaic, ...
            theMidgetRGCMosaicResponses, ...
            spatialFrequenciesTested, ...
            orientationsTested);
        pause
    end


end

function analyzePosition(figNo, theRTVFobj, theSourceRGCindex, theCenterConeTypes, theUnfrozenMidgetRGCMosaic, ...
    theMidgetRGCMosaicResponses, spatialFrequenciesTested, orientationsTested)

        theLcenterVisualRFmap = theRTVFobj.visualRFfromRetinalConePooling(...
                theRTVFobj.LconeRFcomputeStruct.modelConstants, ...
                theRTVFobj.LconeRFcomputeStruct.retinalConePoolingParams.finalValues);
        theLcenterSTFdata = theRTVFobj.visualRFmapPropertiesFromCronerKaplanAnalysis(theLcenterVisualRFmap);
    
        theMcenterVisualRFmap = theRTVFobj.visualRFfromRetinalConePooling(...
                theRTVFobj.MconeRFcomputeStruct.modelConstants, ...
                theRTVFobj.MconeRFcomputeStruct.retinalConePoolingParams.finalValues);
        theMcenterSTFdata = theRTVFobj.visualRFmapPropertiesFromCronerKaplanAnalysis(theMcenterVisualRFmap);


        theSourceRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theSourceRGCindex));
        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = highestExtensionSTF(theSourceRGCResponses, spatialFrequenciesTested, orientationsTested);
    
    
        hFig=figure(figNo); clf;
        set(hFig, 'Position', [10 10 700, 700]);
        
        p1 = plot(spatialFrequenciesTested, allMeasuredSTFs, 'k-'); hold on
        if (theCenterConeTypes == cMosaic.LCONE_ID)
            p2 = plot(spatialFrequenciesTested, theMeasuredSTF, 'ro', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);
        else
            p2 = plot(spatialFrequenciesTested, theMeasuredSTF, 'go', 'MarkerFaceColor', [0.5 1 0.5], 'LineWidth', 1.5);
        end

        p3 = plot(theLcenterSTFdata.spatialFrequencySupport, theLcenterSTFdata.visualSTF/max(theLcenterSTFdata.visualSTF), 'r-', 'LineWidth', 1.5);
        p4 = plot(theMcenterSTFdata.spatialFrequencySupport, theMcenterSTFdata.visualSTF/max(theMcenterSTFdata.visualSTF), 'g--', 'LineWidth', 1.5);
        legend([p1(1) p2 p3 p4], {'measured STFs (optics at (0,0)', 'highest extension measured STF', 'RTVF model STF (L-center)', 'RTVF model STF (M-center)'});

        set(gca, 'XScale', 'log');
        title(sprintf('source RGCpos: %2.3f, %2.3f', ...
             theUnfrozenMidgetRGCMosaic.rgcRFpositionsDegs(theSourceRGCindex,1), ...
             theUnfrozenMidgetRGCMosaic.rgcRFpositionsDegs(theSourceRGCindex,2)));
        theCenterConeTypes
        set(gca, 'FontSize', 16);
        drawnow;

end

function [theSourceRGCindex, theCenterConeTypes] = ...
    sourceRGCindexForRTVFobj(theRTVF, rgcRFcenterConeConnectivityMatrix, coneTypes)
    
    theModelCenterConeIndices = theRTVF.LconeRFcomputeStruct.modelConstants.indicesOfCenterCones;
    for iCone = 1:numel(theModelCenterConeIndices)
        theConeIndex = theModelCenterConeIndices(iCone);
        weights = full(squeeze(rgcRFcenterConeConnectivityMatrix(theConeIndex,:)));
        rgcIndices(iCone) = find(weights > 0);
        theCenterConeTypes(iCone) = coneTypes(theConeIndex);
    end
    theSourceRGCindex = unique(rgcIndices);
    
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