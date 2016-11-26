%% Show the impulse response
%
% The impulse is a flash on a steady background.  You can set the steady
% background and the duration of and such here.
%
% 

%%
ieInit;
clear variables

%% Impulse on a 1 ms time axis.
integrationTime = 5e-3;
tSamples = 100;

sampleTimes = (1:tSamples)*integrationTime;
weights = zeros(tSamples,1);
weights(20:21) = 1;

% Scene parameters in general
sparams.fov       = 0.5;   % Half a degree
sparams.luminance = 50;    % Uniform scene luminance (cd/m2)

% Creates the impulse.  Steady background of 50 cd/m2, then a flash at 100
% cd/m2 for 5 ms, then back to 50.
oiImpulse = oisCreate('impulse','add',weights,...
    'sparams',sparams,...
    'sampleTimes',sampleTimes);

% oiImpulse.visualize;

%% Set the cone mosaic parameters

cMosaic = coneMosaic;    % Default is osLinear
cMosaic.integrationTime = integrationTime;
cMosaic.setSizeToFOV(oiGet(oiImpulse.oiFixed,'fov')*0.8);
cMosaic.emGenSequence(tSamples);

%%  Compute the absorptions

cMosaic.compute(oiImpulse);
cMosaic.window;

%% Compute the current and get the interpolated filters

% The difference between noise flag true and false is actually interesting.
% When the noise is true, the individual cones have relatively little
% reliability.
cMosaic.os.noiseFlag = false;
interpFilters = cMosaic.computeCurrent;
cMosaic.window;

%% If you want to see the interpolated and complete impulse response ..

cMosaic.plot('os impulse response');
hold on;
plot(cMosaic.timeAxis, interpFilters,'o');
grid on; xlabel('Time (sec)');

%%
