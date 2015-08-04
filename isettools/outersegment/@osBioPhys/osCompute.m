function osCompute(obj, sensor, varargin)
    
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());

    if size(varargin{1,1})==0
        params = cell(1,1);
    else
         params = (varargin{1,1});
    end
    
    p = riekeInit;
    expTime = sensorGet(sensor,'exposure time');
    sz = sensorGet(sensor,'size');
    
    % absRate = sensorGet(sensor,'absorptions per second');
    pRate = sensorGet(sensor, 'photon rate');
    
    % Compute background adaptation parameters
    if isfield(params{1,1},'bgVolts')
        bgVolts = params{1,1}.bgVolts; % need to make this an input parameter!
    else
        bgVolts = 0;
    end
    bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
    
    initialState = riekeAdaptSteadyState(bgR, p, sz);
    
    initialState.timeInterval = sensorGet(sensor, 'time interval');
    obj.ConeCurrentSignal  = riekeAdaptTemporal(pRate, initialState);
    
    if isfield(params{1,1},'dc')
        nSamples = length(obj.ConeCurrentSignal);
        obj.ConeCurrentSignal = obj.ConeCurrentSignal - obj.ConeCurrentSignal(:, :, nSamples);
    end
    
    if obj.noiseflag == 1
%         coneSamplingRate = 825;
        paramsNoise.sampTime = 1/sensorGet(sensor, 'time interval');
        obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal*paramsNoise.sampTime)./paramsNoise.sampTime;

%         obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal, paramsNoise);        
        close;
        
    end
end

