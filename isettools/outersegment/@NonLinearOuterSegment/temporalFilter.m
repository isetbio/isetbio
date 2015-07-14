function temporalFilter(obj, sensor, param, varargin)
    
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
%     %
%     % call the parent class method to compute the filter
%     obj.computeFilter();
%     
%     % compute a linear filtering operation
%     obj.outputSignal = conv(obj.filterKernel, obj.inputSignal);
%     
%     % half-wave rectification
%     obj.outputSignal(obj.outputSignal<0) = 0;
%     
%     % exponentiation
%     obj.outputSignal = obj.outputSignal .^ obj.exponent;
    
    
    p = riekeInit;
    expTime = sensorGet(sensor,'exposure time');
    sz = sensorGet(sensor,'size');
    
    % absRate = sensorGet(sensor,'absorptions per second');
    pRate = sensorGet(sensor, 'photon rate');
    
    % Compute background adaptation parameters
    bgVolts = 1; % need to make this an input parameter!
    bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
    
    initialState = riekeAdaptSteadyState(bgR, p, sz);
    obj.ConeCurrentSignal  = riekeAdaptTemporal(pRate, initialState);
    
    if obj.noiseflag == 1
        obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal);
    end
end

