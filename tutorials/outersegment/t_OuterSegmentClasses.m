function t_OuterSegmentClasses

    % generate input signal
    params = paramsGaborColorOpponent();
    [scene, display] = sceneHorwitzHass(params); 
    oi  = oiCreate('wvf human');
    sensor = sensorHorwitzHass(params, scene, oi, display);
    displayClose;
    
    % instantiate a LinearOuterSegment class
    linearOS = osLinear('noiseflag', 1); %osLinear
    
    % instantiate a NonLinearOuterSegment class
    nonlinearOS = osBioPhys('noiseflag', 1); % osBioPhys
    
    % print methods and properties (when verbosity > 0)
    verbosity = 0;
    printClassInfo(verbosity, linearOS);
    printClassInfo(verbosity, nonlinearOS);
    
    % compute linear outer segment response
    linearOS.compute(sensor,params);
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
    