function compute(obj, ...
    initialGridRetinalConePoolingParamsStruct, ...
    multifocalRTVFindicesToBeComputed, ...
    computeLconeCenterComputeStruct, ...
    computeMconeCenterComputeStruct, ...
    exportsDirectory)


    % Allocate memory
    multifocalRTVFobjectsNum = numel(obj.conesNumPooledByTheRFcenterGrid);
    obj.RTVFTobjList = cell(multifocalRTVFobjectsNum ,1);

    if (ischar(multifocalRTVFindicesToBeComputed)) && (strcmp(multifocalRTVFindicesToBeComputed, 'all'))
        multifocalRTVFindicesToBeComputed = 1:multifocalRTVFobjectsNum;
    end


    for iMultifocalRTVFobjIndex = 1:multifocalRTVFobjectsNum

        opticalPositionDegs = obj.opticalPositionGrid(iMultifocalRTVFobjIndex,:);
        conesNumPooled = obj.conesNumPooledByTheRFcenterGrid(iMultifocalRTVFobjIndex);

        if (~ismember(iMultifocalRTVFobjIndex, multifocalRTVFindicesToBeComputed))
             % Report what we are skipping
            fprintf(2,'\n============================================================= \n');
            fprintf(2,'SKIPPING fitting at eccentricity position (degs): (%2.2f, %2.2f) for %d cone RF centers [fit index: %d of %d]', ...
                opticalPositionDegs(1), opticalPositionDegs(2), conesNumPooled, iMultifocalRTVFobjIndex, multifocalRTVFobjectsNum);
            fprintf(2,'\n============================================================= \n');

            continue;
        end
            
        % Report what we are about to do
        fprintf(2,'\n============================================================= \n');
        fprintf(2,'Fitting at eccentricity position (degs): (%2.2f, %2.2f) for a %d cone RF center [fit index: %d of %d]', ...
                opticalPositionDegs(1), opticalPositionDegs(2), conesNumPooled, iMultifocalRTVFobjIndex, multifocalRTVFobjectsNum);
        fprintf(2,'\n============================================================= \n');

        
        % The full Optics params struct
        temporalEquivalentEccDegs = obj.theRGCMosaic.temporalEquivalentEccentricityForEccentricity(opticalPositionDegs);
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        [psfUpsampleFactor, wavefrontSpatialSamples] = psfAndWavefrontSamples(radialTemporalEquivalentEccDegs);

        theOpticsParams = struct(...
            'ZernikeDataBase', obj.opticsParams.ZernikeDataBase, ...
            'examinedSubjectRankOrder', obj.opticsParams.subjectRankOrder, ...
            'pupilDiameterMM', obj.opticsParams.pupilDiameterMM, ...
            'positionDegs', opticalPositionDegs, ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', obj.theRGCMosaic.whichEye, ...
            'subjectRankingEye', 'right eye', ...
            'psfUpsampleFactor', psfUpsampleFactor, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples ...
        );


        % The full target visual RF params
        % Retrieve indices and weights of cones to the targetRGCindex for this RTVF
        theTargetRGCindex = obj.targetRGCindex(iMultifocalRTVFobjIndex);
        indicesOfConesPooledByTheRFcenter = find(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex)> 0);
        weightsOfConesPooledByTheRFcenter = reshape(...
            full(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(indicesOfConesPooledByTheRFcenter,theTargetRGCindex)), ...
            [1 numel(indicesOfConesPooledByTheRFcenter)]);
        
        theTargetVisualRFDoGparams = struct(...
            'retinalConePoolingModel', obj.rfModelParams.retinalConePoolingModel, ...
            'centerConnectableConeTypes', obj.rfModelParams.centerConnectableConeTypes, ...
            'surroundConnectableConeTypes', obj.rfModelParams.surroundConnectableConeTypes, ...
            'indicesOfConesPooledByTheRFcenter', indicesOfConesPooledByTheRFcenter, ...
            'weightsOfConesPooledByTheRFcenter', weightsOfConesPooledByTheRFcenter, ...
            'surroundToCenterRcRatio', obj.visualSTFSurroundToCenterRcRatioGrid(iMultifocalRTVFobjIndex), ...
            'surroundToCenterIntegratedSensitivityRatio', obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid(iMultifocalRTVFobjIndex));


        if (isempty(obj.multiStartsNumRetinalPooling))
            if (conesNumPooled < 3)
                multiStartsNumRetinalPooling = 8*2;
            else
                multiStartsNumRetinalPooling = 4*2;
            end
        else
            multiStartsNumRetinalPooling = obj.multiStartsNumRetinalPooling;
        end

        
        % Check if there are initialRetinalConePoolingParams for this RTVFobj
        % (i.e. at this spatial position, and with this many center cones)
        initialRetinalConePoolingParamsForThisRTVFobj = [];
        
        if (~isempty(initialGridRetinalConePoolingParamsStruct))
            % See if a previous initial retinal cone pooling params exists
            theDictionaryKeys = keys(initialGridRetinalConePoolingParamsStruct.dictionary);
            theTargetKey = sprintf('%d center cones at (%2.3f,%2.3f)', ...
                conesNumPooled, opticalPositionDegs(1), opticalPositionDegs(2));

            if (ismember(theTargetKey, theDictionaryKeys))
                initialRetinalConePoolingParamsForThisRTVFobj = initialGridRetinalConePoolingParamsStruct.dictionary(theTargetKey);
                fprintf(2,'Using previous initial retinal cone pooling params (key: %s)', theTargetKey);
            else
                fprintf(2,'Did not find previous initial retinal cone pooling params for current position. Will use standard initial params\n');
            end
        end

        % Go !
        obj.RTVFTobjList{iMultifocalRTVFobjIndex} = RTVF(...
                obj.theRGCMosaic.inputConeMosaic, ...
                theOpticsParams, ...
                theTargetVisualRFDoGparams, ...
                'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
                'initialRetinalConePoolingParamsStruct', initialRetinalConePoolingParamsForThisRTVFobj, ...
                'computeLconeCenterComputeStruct', computeLconeCenterComputeStruct, ...
                'computeMconeCenterComputeStruct', computeMconeCenterComputeStruct, ...
                'visualizeSpectrallyWeightedPSFs', ~true, ...
                'computedRTVFobjectExportDirectory', exportsDirectory);


        

        theRTFVTobjList = obj.RTVFTobjList;
        theOpticsPositionGrid = obj.opticalPositionGrid;
        theConesNumPooledByTheRFcenterGrid = obj.conesNumPooledByTheRFcenterGrid;
        theVisualSTFSurroundToCenterRcRatioGrid = obj.visualSTFSurroundToCenterRcRatioGrid;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = obj.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

        % Append the grid info to the file saved by the obj.RTVFTobjList{iMultifocalRTVFobjIndex}
        theCurrentRTVFobj = obj.RTVFTobjList{iMultifocalRTVFobjIndex};
        save(theCurrentRTVFobj.computedObjDataFileName, ...
            'iMultifocalRTVFobjIndex', ...
            'theConesNumPooledByTheRFcenterGrid', ...
            'theOpticsPositionGrid', ...
            'theVisualSTFSurroundToCenterRcRatioGrid', ...
            'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
            '-append');

        % Save the currently computed RTVFist to the temporary file
        R2VFTobjFileName = 'tmpRTVFobjects.mat';
        save(R2VFTobjFileName, 'theRTFVTobjList', ...
                        'theOpticsPositionGrid', ...
                        'theConesNumPooledByTheRFcenterGrid', ...
                        'theVisualSTFSurroundToCenterRcRatioGrid', ...
                        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                        '-v7.3');
        fprintf(2,'\n------> Saved current theRTFVTobjList to %s\n', R2VFTobjFileName );
        clear 'theRTFVTobjList';

    end % for iMultifocalRTVFobjIndex 


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