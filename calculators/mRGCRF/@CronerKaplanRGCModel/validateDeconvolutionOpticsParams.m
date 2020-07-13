function validateDeconvolutionOpticsParams(obj,deconvolutionOpticsParams)

    assert(isstruct(deconvolutionOpticsParams), ...
        'deconvolutionOpticsParams must be a struct');
    
    fNames = fieldnames(deconvolutionOpticsParams);
    for k = 1:numel(fNames)
        assert(ismember(fNames{k}, {'PolansWavefrontAberrationSubjectIDsToAverage', 'PolansWavefrontAberrationSubjectIDsToCompute', ...
                                    'quadrantsToAverage', 'quadrantsToCompute'}), ...
           sprintf('field ''%s'' is deconvolutionOpticsParams struct is not valid.', fNames{k}));
    end
   
    if (isfield(deconvolutionOpticsParams, 'quadrantsToAverage'))
        for k = 1:numel(deconvolutionOpticsParams.quadrantsToAverage)
            quadrantName = deconvolutionOpticsParams.quadrantsToAverage{k};
            assert(ismember(quadrantName, obj.validPolansQuadrants), ...
                sprintf('Quadrant ''%s'' is not a valid quadrant name for Polans Optics', deconvolutionOpticsParams.quadrantsToAverage{k}));
        end
    else
        for k = 1:numel(deconvolutionOpticsParams.quadrantsToCompute)
            quadrantName = deconvolutionOpticsParams.quadrantsToCompute{k};
            assert(ismember(quadrantName, obj.validPolansQuadrants), ...
                sprintf('Quadrant ''%s'' is not a valid quadrant name for Polans Optics', deconvolutionOpticsParams.quadrantsToCompute{k}));
         end
    end
    
    if (isfield(deconvolutionOpticsParams, 'PolansWavefrontAberrationSubjectIDsToAverage'))
        for k = 1:numel(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage)
            assert(ismember(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage(k), obj.validPolansSubjectIDs), ...
                sprintf('Subject %d is not a valid subject ID for Polans Optics', deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage(k)));
        end
    else
        for k = 1:numel(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute)
            assert(ismember(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(k), obj.validPolansSubjectIDs), ...
                sprintf('Subject %d is not a valid subject ID for Polans Optics', deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(k)));
        end
    end
    
end

