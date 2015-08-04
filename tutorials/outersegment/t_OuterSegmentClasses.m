function t_OuterSegmentClasses
    % % generate input signal
    % params = paramsGaborColorOpponent();
    % [scene, display] = sceneHorwitzHass(params);
    % oi  = oiCreate('wvf human');
    % sensor = sensorHorwitzHass(params, scene, oi, display);
    % displayClose;
    
    % generate input signal
    %  set up parameters for stimulus
    nSamples = 2000;        % 2000 samples
    timeStep = 1e-4;        % time step
    flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)
    
    %  create human sensor
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
    sensor = sensorSet(sensor, 'time interval', timeStep);
    
    %  create stimulus
    stimulus = zeros(nSamples, 1);
    stimulus(1) = flashIntens;
    stimulus = reshape(stimulus, [1 1 nSamples]);
    
    % set photon rates
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % compute model current and baseline correct
    params.bgVolts  = 0;
    % params = 1;
    
    % instantiate a LinearOuterSegment class
    linearOS = osLinear('noiseflag', 1); %osLinear
    
    % instantiate a NonLinearOuterSegment class
    nonlinearOS = osBioPhys('noiseflag', 1); % osBioPhys
    
    % print methods and properties (when verbosity > 0)
    verbosity = 0;
    printClassInfo(verbosity, linearOS);
    printClassInfo(verbosity, nonlinearOS);
    
    % compute linear outer segment response
    linearOS.compute(sensor);
    linearOS.get('noiseflag')
    
    % compute nonlinear outer segment response
    nonlinearOS.compute(sensor,params);
    linearOS.get('noiseflag')
    
    % plot results
    linearOS.plotResults(sensor);
    nonlinearOS.plotResults(sensor);
    
end


function printClassInfo(verbosity, os)
    if (verbosity > 0)
        fprintf('<strong>\n\nHit enter to list the public methods of %s</strong>', class(os));
        pause;
        methods(os)

        fprintf('<strong>\n\nHit enter to list the public properties of %s</strong>', class(os));
        pause; 
        properties(os)
    end
end
    