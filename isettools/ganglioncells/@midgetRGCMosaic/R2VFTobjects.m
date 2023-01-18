function [RTVFTobjList, ...
          theSamplingPositionGrid, ...
          theConesNumPooledByTheRFcenterGrid, ...
          theVisualSTFSurroundToCenterRcRatioGrid, ...
          theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = R2VFTobjects(...
                    RTVobjIndicesToBeComputed, ...
                    theMidgetRGCMosaic, ...
                    eccentricitySamplingGrid, ...
                    mosaicSurroundParams, opticsParams, fitParams)

    ZernikeDataBase = opticsParams.ZernikeDataBase;
    subjectRankOrder = opticsParams.subjectRankOrder;
    pupilDiameterMM = opticsParams.pupilDiameterMM; 

    centerConnectableConeTypes = mosaicSurroundParams.centerConnectableConeTypes;
    surroundConnectableConeTypes = mosaicSurroundParams.surroundConnectableConeTypes;
    coneWeightsCompensateForVariationsInConeEfficiency = mosaicSurroundParams.coneWeightsCompensateForVariationsInConeEfficiency;
    visualRFmodel = mosaicSurroundParams.visualRFmodel;
    retinalConePoolingModel = mosaicSurroundParams.retinalConePoolingModel;
    targetSTFmatchMode = mosaicSurroundParams.targetSTFmatchMode;

    multiStartsNumDoGFit = fitParams.multiStartsNumDoGFit;
    exportsDirectory = fitParams.exportsDirectory;

    % Find out the range of cones in the RF center
    allConesNumPooledByTheRFcenters = full(sum(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix,1));
    conesNumPooledByTheRFcenters = unique(allConesNumPooledByTheRFcenters);
    fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenters);

    % Generate grid of optical positions x conesNumInRFcenter
    [RTVTOpticalPositionGrid, RTVTConesNumInRFcenterGrid] = ...
        meshgrid(1:size(eccentricitySamplingGrid,1), conesNumPooledByTheRFcenters);
    
    % Allocate memory
    RTVTobjectsNum = numel(RTVTConesNumInRFcenterGrid);
    RTVFTobjList = cell(RTVTobjectsNum,1);
    theSamplingPositionGrid = zeros(RTVTobjectsNum,2);
    theConesNumPooledByTheRFcenterGrid = zeros(RTVTobjectsNum,1);
    theVisualSTFSurroundToCenterRcRatioGrid = zeros(RTVTobjectsNum,1);
    theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = zeros(RTVTobjectsNum,1);

    if (ischar(RTVobjIndicesToBeComputed)) && (strcmp(RTVobjIndicesToBeComputed, 'all'))
        RTVobjIndicesToBeComputed = 1:RTVTobjectsNum;
    end

    for iRTVobjIndex = 1:RTVTobjectsNum
        % Sampling position (within the  mosaic)
        eccPositionIndex = RTVTOpticalPositionGrid(iRTVobjIndex);
        samplingPositionDegs = eccentricitySamplingGrid(eccPositionIndex,:);

        % Cones num in RF center
        conesNumPooled = RTVTConesNumInRFcenterGrid(iRTVobjIndex);

        if (conesNumPooled < 3)
            multiStartsNum = 8;
        else
            multiStartsNum = 4;
        end

        % Get indices of center cones for the RGC that is closest to the samplingPositionDegs
        indicesOfRGCsWithThisManyCenterCones = find(allConesNumPooledByTheRFcenters == conesNumPooled);
        [~,idx] = min(sum((bsxfun(@minus, theMidgetRGCMosaic.rgcRFpositionsDegs(indicesOfRGCsWithThisManyCenterCones,:), samplingPositionDegs)).^2,2));
        theTargetRGCindex = indicesOfRGCsWithThisManyCenterCones(idx);

        % Update the samplingPositionDegs to reflect the actual position of theTargetRGCindex
        samplingPositionDegs = theMidgetRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);

        % Check for initialRetinalConePoolingParams
        initialRetinalConePoolingParamsStruct = [];
        if (isfield(fitParams, 'initialRetinalConePoolingParamsStruct')) && ...
           (~isempty(fitParams.initialRetinalConePoolingParamsStruct))
            % See if a previous initial retinal cone pooling params exists
            theDictionaryKeys = keys(fitParams.initialRetinalConePoolingParamsStruct.dictionary);
            theTargetKey = sprintf('%d center cones at (%2.3f,%2.3f)', ...
                conesNumPooled, samplingPositionDegs(1), samplingPositionDegs(2));

            if (ismember(theTargetKey, theDictionaryKeys))
                initialRetinalConePoolingParamsStruct = fitParams.initialRetinalConePoolingParamsStruct.dictionary(theTargetKey);
            else
                fprintf(2,'Did not find previous initial retinal cone pooling params for current position. Will use standard initial params\n');
            end
        end

        indicesOfConesPooledByTheRFcenter = find(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex)> 0);
        typesOfConesPooledByTheRFcenter = theMidgetRGCMosaic.inputConeMosaic.coneTypes(indicesOfConesPooledByTheRFcenter);

        % Assert that we have the correct number of center cones
        assert((numel(indicesOfConesPooledByTheRFcenter) == conesNumPooled), ...
            sprintf('indicesOfConesPooledByTheRFcenter should have %d entries but it has %d.', numel(indicesOfConesPooledByTheRFcenter), conesNumPooled));
        
        % Assert that the center cones are all connectable
        assert(all(ismember(typesOfConesPooledByTheRFcenter, centerConnectableConeTypes)), ...
            sprintf('indicesOfConesPooledByTheRFcenter are not all connectable'));

        % Unit weights for all center cones
        weightsOfConesPooledByTheRFcenter = ones(1,numel(indicesOfConesPooledByTheRFcenter));

        % Rs/Rc ratio
        surroundToCenterRcRatio = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

        % Temporal-equivalent eccentricity based SCint sensitivity ratio
        temporalEquivalentEccDegs = theMidgetRGCMosaic.temporalEquivalentEccentricityForEccentricity(samplingPositionDegs);
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        scIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);
        
        % Optics params
        [psfUpsampleFactor, wavefrontSpatialSamples] = psfAndWavefrontSamples(radialTemporalEquivalentEccDegs);
        theOpticsParams = struct(...
            'ZernikeDataBase', ZernikeDataBase, ...
            'examinedSubjectRankOrder', subjectRankOrder, ...
            'pupilDiameterMM', pupilDiameterMM, ...
            'positionDegs', samplingPositionDegs, ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', theMidgetRGCMosaic.whichEye, ...
            'subjectRankingEye', 'right eye', ...
            'psfUpsampleFactor', psfUpsampleFactor, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples ...
        );

      
        theTargetVisualRFDoGparams = struct(...
            'visualRFmodel', visualRFmodel, ... 
            'retinalConePoolingModel', retinalConePoolingModel, ...
            'centerConnectableConeTypes', centerConnectableConeTypes, ...
            'surroundConnectableConeTypes', surroundConnectableConeTypes, ...
            'coneWeightsCompensateForVariationsInConeEfficiency', coneWeightsCompensateForVariationsInConeEfficiency,  ...
            'indicesOfConesPooledByTheRFcenter', indicesOfConesPooledByTheRFcenter, ...
            'weightsOfConesPooledByTheRFcenter', weightsOfConesPooledByTheRFcenter, ...
            'conesNumPooledByTheRFcenter', conesNumPooled, ...
            'surroundToCenterRcRatio', surroundToCenterRcRatio, ...
            'surroundToCenterIntegratedSensitivityRatio', scIntSensitivity);

        

        theSamplingPositionGrid(iRTVobjIndex,:) = samplingPositionDegs;
        theConesNumPooledByTheRFcenterGrid(iRTVobjIndex) = conesNumPooled;
        theVisualSTFSurroundToCenterRcRatioGrid(iRTVobjIndex) = surroundToCenterRcRatio;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid(iRTVobjIndex) = scIntSensitivity;

        if (ismember(iRTVobjIndex, RTVobjIndicesToBeComputed))

            % Report what we are about to do
            fprintf(2,'\n============================================================= \n');
            fprintf(2,'Fitting at eccentricity position (degs): (%2.2f, %2.2f) for a %d cone RF center [fit index: %d of %d]', ...
                samplingPositionDegs(1), samplingPositionDegs(2), numel(typesOfConesPooledByTheRFcenter), iRTVobjIndex, RTVTobjectsNum);
            fprintf(2,'\n============================================================= \n');

            % Compute the RetinaToVisualFieldTransformer for this (optics & cone position, conesNumPooled) pair
            % We are simulating the Croner&Kaplan estimation (1D STF, and we
            % are matching the fitted STF DoG model params)
            RTVFTobjList{iRTVobjIndex} = RetinaToVisualFieldTransformer(...
                theMidgetRGCMosaic.inputConeMosaic, ...
                theOpticsParams, ...
                theTargetVisualRFDoGparams, ...
                'simulateCronerKaplanEstimation', true, ...
                'targetSTFmatchMode', targetSTFmatchMode, ...
                'initialRetinalConePoolingParamsStruct', initialRetinalConePoolingParamsStruct, ...
                'multiStartsNumDoGFit', multiStartsNumDoGFit, ...
                'multiStartsNum', multiStartsNum , ...
                'doDryRunFirst', true, ...
                'computedRTVObjectExportDirectory',  exportsDirectory);

            % Append the grids and the current iRTVobjIndex
            save(RTVFTobjList{iRTVobjIndex}.computedObjDataFileName, ...
                'iRTVobjIndex', ...
                'theSamplingPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                '-append');
        else
            % Report what we are skipping
            fprintf(2,'\n============================================================= \n');
            fprintf(2,'SKIPPING fitting at eccentricity position (degs): (%2.2f, %2.2f) for a %d cone RF center [fit index: %d of %d]', ...
                samplingPositionDegs(1), samplingPositionDegs(2), numel(typesOfConesPooledByTheRFcenter), iRTVobjIndex, RTVTobjectsNum);
            fprintf(2,'\n============================================================= \n');

          
        
        end
    end
end

function [psfUpsampleFactor, wavefrontSpatialSamples] = psfAndWavefrontSamples(radialEccDegs)
    if (radialEccDegs <= 0.5)
        % Cones are tiny, so upsample the spatial resolution
        psfUpsampleFactor = 2;
        wavefrontSpatialSamples = 301;
    elseif (radialEccDegs <= 8)
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 401;
    elseif (radialEccDegs <= 14)
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 601;   
    else
        psfUpsampleFactor = 1;
        wavefrontSpatialSamples = 701;
    end
end
