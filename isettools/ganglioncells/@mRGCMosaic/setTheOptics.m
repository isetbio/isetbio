function  setTheOptics(obj, opticsParams)

    % Generate the native optics.
    % These optics form the basis on which surround cone weights are optimized
    % for wiring so as to generate RF with the target visual properties
    dataOut = obj.generateOptics(opticsParams);

    if (isempty(opticsParams.positionDegs))
        % Save the native optics params and the native optics
        obj.theNativeOpticsParams = dataOut.opticsParams;
        obj.theNativeOptics = dataOut.theOptics;
        fprintf('Generated native optics at mosaic''s center\n');
    else
        obj.theCustomOpticsParams = dataOut.opticsParams;
        obj.theCustomOptics = dataOut.theOptics;
        fprintf('Generated custom optics at (%2.1f, %2.1f)\n', ...
            opticsParams.positionDegs(1), opticsParams.positionDegs(2));
    end

end


    