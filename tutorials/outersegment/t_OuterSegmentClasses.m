function t_OuterSegmentClasses
% A tutorial function for the outer segment object. The outer segment
% object converts cone isomerizations (R*) over time to cone outersegment 
% current (pA) over time. outerSegment is a parent class for the osLinear 
% and osBioPhys objects. The osLinear object calculates the outer segment 
% current by convolving linear filters for the L, M and S cones with the 
% isomerization signal, while the osBioPhys objects calculates the outer
% segment current with a nonlinear difference equation model by Fred Rieke. 
% This tutorial function illustrates how to initialize the outer segment 
% objects, compute the current signal, add noise to the current signal 
% (according to another model developed by Fred Rieke), and plot results.
%
% Additionally, a new version of coneAdapt is implemented by the
% outerSegment class that allows for backward comptability with the
% previous coneAdapt function. A call of coneAdapt will now return the same
% value as the previous version, but also returns an osLinear or osBioPhys
% object as well. The new coneAdapt stores properties and calls methods
% using the outerSegment object methods.
% 
% In order to have the names of the public properties and methods of each
% subclass printed in the command window, set 'verbosity' equal to 1 below.
% 
% t_OuterSegmentClasses
% 
% 7/2015 JRG NC DHB
    
    % Generate input signal.
    % Set up parameters for stimulus
    nSamples = 2000;        % 2000 samples
    timeStep = 1e-4;        % time step
    flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)
    
    % Create human sensor.
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
    sensor = sensorSet(sensor, 'time interval', timeStep);
    
    % Create stimulus - impulse.
    stimulus = zeros(nSamples, 1);
    stimulus(1) = flashIntens;
    stimulus = reshape(stimulus, [1 1 nSamples]);
    
    % Set photon rates.
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % Instantiate an osLinear class.
    linearOS = osLinear('noiseFlag', 1); %osLinear
    
    % Instantiate an osBioPhys class.
    adaptedOS = osBioPhys(); % osBioPhys
    adaptedOS = osBioPhysSet(adaptedOS, 'noiseFlag', 1);
    
    % Print methods and properties (when verbosity > 0).
    verbosity = 0;
    printClassInfo(verbosity, linearOS);
    printClassInfo(verbosity, adaptedOS);
    
    % Compute linear outer segment response.
    linearOS = osLinearCompute(linearOS, sensor);
    % params.offset = 0;
    % linearOS = osLinearCompute(linearOS, sensor, params);
    osLinearGet(linearOS, 'noiseFlag');
    
    % Compute nonlinear outer segment response.
    adaptedOS = osBioPhysCompute(adaptedOS, sensor);
    % params.bgVolts = 0; params.offset = 0;
    % adaptedOS = osBioPhysCompute(adaptedOS, sensor, params);
    osBioPhysGet(adaptedOS, 'noiseFlag');
    
    % Plot results.
    osLinearPlot(linearOS, sensor);
    osBioPhysPlot(adaptedOS, sensor);
    
    % Alternate function call, provides backwards compatibility with
    % coneAdapt but calls through outerSegment object code.
    % osCurrent = coneAdaptAlt(sensor, 'linear');
    % [osCurrent linearOSalt] = coneAdaptAlt(sensor, 'linear');
    % osCurrent = coneAdaptAlt(sensor, 'rieke');
    % [osCurrent adaptedOSalt] = coneAdaptAlt(sensor, 'rieke');
    
end


function printClassInfo(verbosity, os)
    if (verbosity > 0)
        fprintf('<strong>\n\nPublic methods of %s:</strong>', class(os));
        methods(os)

        fprintf('<strong>\n\nPublic properties of %s:</strong>', class(os));
        properties(os)
    end
end
    