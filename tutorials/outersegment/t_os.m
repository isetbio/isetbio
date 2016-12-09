% t_os
%
% A tutorial function for the outer segment object. The outer segment
% object converts cone isomerizations (R*) over time to cone outersegment
% current (pA) over time. outerSegment is a parent class for the osLinear
% and osBioPhys objects. 
% 
% The osLinear object calculates the outer segment
% current by convolving linear filters for the L, M and S cones with the
% isomerization signal, while the osBioPhys objects calculates the outer
% segment current with a nonlinear difference equation model by Fred Rieke.
% This tutorial function illustrates how to initialize the outer segment
% objects, compute the current signal, add noise to the current signal
% (according to another model developed by Fred Rieke), and plot results.
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
sensor = sensorSet(sensor, 'photon rate', stimulus);

%% Instantiate an osLinear object.
linearOS = osCreate('linear');%,'noiseFlag', 1);
linearOS = osSet(linearOS, 'noise flag', 'none');
% linearOS = osLinear('noiseFlag', 1); %osLinear

%% Compute linear outer segment response.
linearOS = osCompute(linearOS, sensor);
osGet(linearOS, 'noise flag');
linearOS.plot('input');

%% Instantiate an osBioPhys object.
adaptedOS = osCreate('biophys');
% adaptedOS = osBioPhys(); % osBioPhys
adaptedOS = osSet(adaptedOS, 'noiseFlag', 0);

%% Compute nonlinear outer segment response.
adaptedOS = osCompute(adaptedOS, sensor);
osGet(adaptedOS, 'noiseFlag');

%% Plot results.
osPlot(linearOS, sensor);
osPlot(adaptedOS, sensor);

%%


% Set photons.
sensor = sensorSet(sensor, 'photon rate', stimulus);

%% Instantiate an osLinear object.
identityOS = osCreate('identity');%,'noiseFlag', 1);
identityOS = osSet(identityOS, 'noiseFlag', 0);

%% Compute linear outer segment response.
identityOS = osCompute(identityOS, sensor);
osGet(identityOS, 'noiseFlag');

%% Plot results.
osPlot(identityOS, sensor);
