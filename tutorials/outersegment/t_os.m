% t_os
%
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
% In order to have the names of the public properties and methods of each
% subclass printed in the command window, set 'verbosity' equal to 1 below.
%
% 7/2015 JRG NC DHB

%% Generate input signal.
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

% Set photons.
sensor = sensorSet(sensor, 'photons', stimulus);

%% Instantiate an osLinear class.
linearOS = osCreate('linear');%,'noiseFlag', 1);
linearOS = osSet(linearOS,'noiseFlag',1);
% linearOS = osLinear('noiseFlag', 1); %osLinear

%% Instantiate an osBioPhys class.
adaptedOS = osCreate('biophys');
% adaptedOS = osBioPhys(); % osBioPhys
adaptedOS = osSet(adaptedOS, 'noiseFlag', 1);

%% Compute linear outer segment response.
linearOS = osCompute(linearOS, sensor);
% params.offset = 0;
% linearOS = osLinearCompute(linearOS, sensor, params);
osGet(linearOS, 'noiseFlag');

%% Compute nonlinear outer segment response.
adaptedOS = osCompute(adaptedOS, sensor);
% params.bgVolts = 0; params.offset = 0;
% adaptedOS = osBioPhysCompute(adaptedOS, sensor, params);
osGet(adaptedOS, 'noiseFlag');

%% Plot results.
osPlot(linearOS, sensor);
osPlot(adaptedOS, sensor);

%% Print methods and properties (when verbosity > 0).
verbosity = 1;
if (verbosity > 0)
    fprintf('<strong>\n\nPublic methods of %s:</strong>', class(linearOS));
    methods(linearOS)    
    fprintf('<strong>\n\nPublic properties of %s:</strong>', class(linearOS));
    properties(linearOS)    
    fprintf('<strong>\n\nPublic methods of %s:</strong>', class(adaptedOS));
    methods(adaptedOS)    
    fprintf('<strong>\n\nPublic properties of %s:</strong>', class(adaptedOS));
    properties(adaptedOS)
end