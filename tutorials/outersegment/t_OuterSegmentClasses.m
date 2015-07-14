function demoOuterSegmentClasses

    % generate input signal
    params = gabor_color_opponent_params();
    [scene, display] = build_scene_horwitz_hass_2015(params); 
    oi  = oiCreate('wvf human');
    sensor = build_sensor_horwitz_hass_2015(params, scene, oi, display);
    displayClose;
    
    % instantiate a LinearOuterSegment class
    linearOS = LinearOuterSegment('noiseflag', 1);
    
    % instantiate a NonLinearOuterSegment class
    nonlinearOS = NonLinearOuterSegment('noiseflag', 1);
    
    % print methods and properties (when verbosity > 0)
    verbosity = 1;
    printClassInfo(verbosity, linearOS);
    printClassInfo(verbosity, nonlinearOS);
    
    % compute linear outer segment response
    linearOS.temporalFilter(sensor,params);
    linearOS.get('noiseflag')
    
    % compute nonlinear outer segment response
    nonlinearOS.temporalFilter(sensor,params);
    linearOS.get('noiseflag')
    
    % plot results
%     plotResults(1, dt, linearOS, sensor);
%     plotResults(2, dt, nonlinearOS);  
    
end

function plotResults(figNum, dt, os, sensor)

    h = figure(figNum);
    set(h, 'Name', sprintf('Output of %s', class(os)));
    set(h, 'Position', [10+50*figNum 10+50*figNum, 1024 256]);
    subplot(1,3,1);
    plot((0:numel(os.filterKernel)-1)*dt, os.inputSignal, 'k-');
    title('input signal');
    subplot(1,3,2);
    % plot((0:numel(os.filterKernel)-1)*dt, os.filterKernel, 'r-');
    hold on;
%     plot(os.sConeFilter,'b'); plot(os.mConeFilter,'g'); plot(os.lConeFilter,'r');
    title('filter kernel')
    subplot(1,3,3);
    plot((0:numel(os.outputSignal)-1)*dt, os.outputSignal, 'k-');
    title('output signal');
    drawnow;
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
    