function generateNativeOptics(obj, opticsParams)

    % Generate the native optics.
    % These optics form the basis on which surround cone weights are optimized
    % for wiring so as to generate RF with the target visual properties
    dataOut = obj.generateOptics(opticsParams);

    % Save the native optics params and the native optics
    obj.theNativeOpticsParams = dataOut.opticsParams;
    obj.theNativeOptics = dataOut.theOptics;
end


    