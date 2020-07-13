function validateDeconvolutionOpticsParams(obj,deconvolutionOpticsParams)

    assert(isstruct(deconvolutionOpticsParams), ...
        'deconvolutionOpticsParams must be a struct');
    
    fNames = fieldnames(deconvolutionOpticsParams);
    for k = 1:numel(fNames)
        assert(ismember(fNames{k}, {'PolansWavefrontAberrationSubjectIDsToAverage', 'quadrantsToAverage'}), ...
           sprintf('field ''%s'' is deconvolutionOpticsParams struct is not valid.', fNames{k}));
    end
    
    for k = 1:numel(deconvolutionOpticsParams.quadrantsToAverage)
        assert(ismember(deconvolutionOpticsParams.quadrantsToAverage{k}, obj.validPolansQuadrants), ...
            sprintf('Quadrant ''%s'' is not a valid quadrant name for Polans Optics', deconvolutionOpticsParams.quadrantsToAverage{k}));
    end
    for k = 1:numel(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage)
        assert(ismember(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDToAverage(k), obj.alidPolansSubjectIDs), ...
            sprintf('Subject %d is not a valid subject ID for Polans Optics', deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDToAverage(k)));
    end
end

