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
weights(20:25) = 1;

%
sparams.fov       = 0.5;
sparams.luminance = 50;

oiImpulse = oisCreate('impulse','add',weights,'sparams',sparams,'sampleTimes',sampleTimes);
oiImpulse.visualize;

%%
cMosaic = coneMosaic;
cMosaic.integrationTime = integrationTime;
cMosaic.setSizeToFOV(oiGet(oiImpulse.oiFixed,'fov')*0.8);

% The amount of time for the cone measurements should be the amount of time
% for the oiSequence
T = oiImpulse.timeAxis(end);
cMosaic.emGenSequence(tSamples);

%%
cMosaic.compute(oiImpulse);
cMosaic.window;

%% Compute the current and get the interpolated filters

% Return is only meaningful for osLinear model
interpFilters = cMosaic.computeCurrent;

cMosaic.plot('os impulse response');
hold on;
plot(cMosaic.timeAxis, interpFilters,'o');
grid on; xlabel('Time (sec)');

%%
cMosaic.window;

%%