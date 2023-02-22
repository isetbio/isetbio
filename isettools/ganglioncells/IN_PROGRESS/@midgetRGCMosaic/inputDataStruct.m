function theInputDataStruct = inputDataStruct(computeInputDataType)

    if (~ismember(computeInputDataType, midgetRGCMosaic.validComputeInputDataTypes))
        fprintf(2, '\nValid compute input data types for the @midgetRGCMosaic are: ');
        for i = 1:numel(midgetRGCMosaic.validComputeInputDataTypes)
            fprintf(2,'\n\t ''%s'' ', midgetRGCMosaic.validComputeInputDataTypes{i});
        end
        error('\nThe passed compute input data type (''%s'') is not valid.\n', computeInputDataType);
    end

    switch (computeInputDataType)

        case midgetRGCMosaic.CONE_MOSAIC_RESPONSE_COMPUTE_INPUT_DATA_TYPE
            % skeleton struct
            theInputDataStruct = struct(...
                'type', computeInputDataType, ...
                'inputConeMosaicSpatioTemporalActivation', [], ...
                'inputConeMosaicActivationTemporalSupport', 0 ...
                );

        case midgetRGCMosaic.SCENE_COMPUTE_INPUT_DATA_TYPE
            % skeleton struct
            theInputDataStruct = struct(...
                'type', computeInputDataType, ...
                'theTestScene', struct(), ...    
                'theNullScene', struct(), ...
                'wavefrontOpticsPositionDegs', [0 0], ...
                'opticalImagePositionDegs', 'mosaic-centered', ...
                'normalizeConeResponsesWithRespectToNullScene', false);
            
    end


end