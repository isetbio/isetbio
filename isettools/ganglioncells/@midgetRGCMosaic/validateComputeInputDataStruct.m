function validateComputeInputDataStruct(inputDataStruct)
    assert(isstruct(inputDataStruct), 'midgetRGCMosaic.validateComputeInputDataStruct: inputData must be a struct');
    assert(isfield(inputDataStruct, 'type'), 'midgetRGCMosaic.validateComputeInputDataStruct: inputData struct must have a ''type'' field.');
    
    % Validate computeInputDataType
    if (~ismember(inputDataStruct.type, midgetRGCMosaic.validComputeInputDataTypes))
        fprintf(2, '\nValid compute input data types for the @midgetRGCMosaic are: ');
        for i = 1:numel(midgetRGCMosaic.validComputeInputDataTypes)
            fprintf(2,'\n\t ''%s'' ', midgetRGCMosaic.validComputeInputDataTypes{i});
        end
        error('\nmidgetRGCMosaic.validateComputeInputDataStruct: the passed compute input data type (''%s'') is not valid.\n', computeInputDataType);
    end

    switch (inputDataStruct.type)
        case midgetRGCMosaic.CONE_MOSAIC_RESPONSE_COMPUTE_INPUT_DATA_TYPE
            % Specific validation
            assert(ndims(inputDataStruct.inputConeMosaicSpatioTemporalActivation) == 3, ...
                'midgetRGCMosaic.validateComputeInputDataStruct: inputConeMosaicSpatioTemporalActivation must have 3 dimensions [nTrials x nTimePoints x nCones]');
    
            assert(size(inputDataStruct.inputConeMosaicSpatioTemporalActivation,2) == numel(inputDataStruct.inputConeMosaicActivationTemporalSupport), ...
                'midgetRGCMosaic.validateComputeInputDataStruct: size(coneMosaicResponses,2) (%d) does not equal the length of temporal support (%d)', ...
                size(inputDataStruct.inputConeMosaicSpatioTemporalActivation,2), numel(inputDataStruct.inputConeMosaicActivationTemporalSupport));

        case midgetRGCMosaic.SCENE_COMPUTE_INPUT_DATA_TYPE
            % Specific validation
            assert(isstruct(inputDataStruct.theTestScene), 'midgetRGCMosaic.validateComputeInputDataStruct: inputDataStruct.theTestScene must be a scene');
            assert(isstruct(inputDataStruct.theNullScene), 'midgetRGCMosaic.validateComputeInputDataStruct: inputDataStruct.theNullScene must be a scene');

            assert(...
                 (isnumeric(inputDataStruct.wavefrontOpticsPositionDegs)) && ...
                 (numel(inputDataStruct.wavefrontOpticsPositionDegs) == 2), ...
                 'midgetRGCMosaic.validateComputeInputDataStruct: inputDataStruct.wavefrontOpticsPositionDegs must be a 2-element vector');

            assert(...
                 (ischar(inputDataStruct.opticalImagePositionDegs)) || ...
                 ((isnumeric(inputDataStruct.opticalImagePositionDegs)) && (numel(inputDataStruct.opticalImagePositionDegs) == 2)), ...
                 'midgetRGCMosaic.compute: inputDataStruct.opticalImagePositionDegs must be either a 2-element vector or set to ''mosaic-centered''.');

            assert(islogical(inputDataStruct.normalizeConeResponsesWithRespectToNullScene), 'midgetRGCMosaic.validateComputeInputDataStruct: inputDataStruct.normalizeConeResponsesWithRespectToNullScene must be a boolean');
    end

end
