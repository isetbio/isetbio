function setTheOptics(obj, opticsParams)

    % Generate the native optics.
    % These optics form the basis on which surround cone weights are optimized
    % for wiring so as to generate RF with the target visual properties
    dataOut = obj.generateOptics(opticsParams);

    if (isfield(opticsParams, 'employMonochromaticVlambdaWeightedPSF'))
        if (opticsParams.employMonochromaticVlambdaWeightedPSF)
            dataOut.theOptics = monochromaticVlambdaWeightedOpticsFromOI(dataOut.theOptics);
        end
    end

    if (isempty(opticsParams.positionDegs))
        % Save the native optics params and the native optics
        obj.theNativeOpticsParams = dataOut.opticsParams;
        obj.theNativeOptics = dataOut.theOptics;
        if (opticsParams.examinedSubjectRankOrder == 0)
            fprintf('Generated AO optics at mosaic''s center\n');
        else
            fprintf('Generated native optics at mosaic''s center\n');
        end
    else
        obj.theCustomOpticsParams = dataOut.opticsParams;
        obj.theCustomOptics = dataOut.theOptics;
        fprintf('Generated custom optics at (%2.1f, %2.1f)\n', ...
            opticsParams.positionDegs(1), opticsParams.positionDegs(2));
    end

end    