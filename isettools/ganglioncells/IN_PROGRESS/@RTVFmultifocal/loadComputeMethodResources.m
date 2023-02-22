function loadComputeMethodResources(obj, coneMosaicResponsesFileName)
    % Load cone mosaic STF responses from disk
    obj.stfComputeMethodResources = load(coneMosaicResponsesFileName);
    validateFieldNames(fieldnames(obj.stfComputeMethodResources));

    conesNum = size(obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs,1);
    assert(size(obj.stfComputeMethodResources.theConeMosaicResponses,4) == conesNum, ...
        'size of responses is not consistent with the # of input cones in theRGCMosaic');

    fprintf('\nComputing cone contrast responses ...');
    % Convert cone responses to cone contrast responses
    conesNum = numel(obj.stfComputeMethodResources.coneMosaicNullResponses);
    nullResponses = reshape(obj.stfComputeMethodResources.coneMosaicNullResponses, [1 1 1 conesNum]);

    coneIndicesWithZeroNullResponse = find(nullResponses == 0);
    normalizingResponses = 1./nullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;

    obj.stfComputeMethodResources.theConeMosaicResponses = ...
        bsxfun(@times, bsxfun(@minus, obj.stfComputeMethodResources.theConeMosaicResponses, nullResponses), normalizingResponses);

    % No longer need the null responses
    obj.stfComputeMethodResources = rmfield(obj.stfComputeMethodResources, 'coneMosaicNullResponses');
    fprintf('Done !\n');
end

function validateFieldNames(existingFields)

    requiredFields = {
        'coneContrasts', ...
        'coneMosaicNullResponses', ...
        'opticsPositionDegs', ...
        'orientationsTested', ...
        'spatialFrequenciesTested', ...
        'spatialPhasesDegs', ...
        'theConeMosaicResponses' ...
        }; 

    for iField = 1:numel(existingFields)
        assert(ismember(existingFields{iField}, requiredFields), ...
            sprintf('Fieldname ''%s'', was not expected', existingFields{iField}));
    end

    for iField = 1:numel(requiredFields)
        assert(ismember(requiredFields{iField}, existingFields), ...
            sprintf('Fieldname ''%s'', was not found', requiredFields{iField}));
    end


end


