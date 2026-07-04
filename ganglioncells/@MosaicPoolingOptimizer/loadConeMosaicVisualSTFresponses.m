function loadConeMosaicVisualSTFresponses(obj, responsesFileName)

    load(responsesFileName, 'theNativeOpticsParams', ...
        'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses == 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 1 numel(normalizingResponses)]);
          
    obj.inputConeMosaicVisualSTFdata.orientationsTested = orientationsTested;
    obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested = spatialFrequenciesTested;
    obj.inputConeMosaicVisualSTFdata.spatialPhasesDegs = spatialPhasesDegs;

    % Transform excitations to contrast responses
    theConeMosaicNullResponses = reshape(theConeMosaicNullResponses, [1 1 1 numel(theConeMosaicNullResponses)]);
    obj.inputConeMosaicVisualSTFdata.responseModulations = ...
        bsxfun(@times, bsxfun(@minus, theConeMosaicSTFresponses, theConeMosaicNullResponses), ...
                            normalizingResponses);

    obj.inputConeMosaicVisualSTFdata.metaData.coneContrasts = coneContrasts;
    obj.inputConeMosaicVisualSTFdata.metaData.theNativeOpticsParams = theNativeOpticsParams;

end


