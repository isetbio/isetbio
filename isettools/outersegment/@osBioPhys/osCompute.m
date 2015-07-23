function osCompute(obj, sensor, params, varargin)
    
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());

    p = riekeInit;
    expTime = sensorGet(sensor,'exposure time');
    sz = sensorGet(sensor,'size');
    
    % absRate = sensorGet(sensor,'absorptions per second');
    pRate = sensorGet(sensor, 'photon rate');
    
    % Compute background adaptation parameters
    bgVolts = params.bgVolts; % need to make this an input parameter!
    bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
    
    initialState = riekeAdaptSteadyState(bgR, p, sz);
    
    initialState.timeInterval = sensorGet(sensor, 'time interval');
    obj.ConeCurrentSignal  = riekeAdaptTemporal(pRate, initialState);
    
    if isfield(params,'dc')
        nSamples = length(obj.ConeCurrentSignal);
        obj.ConeCurrentSignal = obj.ConeCurrentSignal - obj.ConeCurrentSignal(:, :, nSamples);
    end
    
    if obj.noiseflag == 1
        obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal);
    end
end

