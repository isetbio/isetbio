%% s_humanPigments
%
% Create a human sensor and plot the spectral QE that combines the lens,
% macular and cone pigments
%
% Plot the lens, macula and cone pigments separatelyk
%
% Adjust the lens and macular densities, then replot the spectral QE
%
% HJ Vistasoft Team, Copyright 2014


%%
s_initISET

%% Show the default human spectral qe

sensor = sensorCreate('human');

vcNewGraphWin;
w = sensorGet(sensor,'wave');
f = sensorGet(sensor,'spectral qe');
plot(w,f);

%%  Show the components - lens transmittance

f = sensorGet(sensor,'human lens transmittance');
plot(w,f);
grid on;
xlabel('Wavelength (nm)')
ylabel('Transmittance')


%% Macular transmittance

f = sensorGet(sensor,'human macular transmittance');
plot(w,f);
grid on;
xlabel('Wavelength (nm)')
ylabel('Transmittance')

%% Set the parameters and plot again

sensor = sensorSet(sensor,'human macular density',0.05);
sensor = sensorSet(sensor,'human lens density',0.5);

vcNewGraphWin;
w = sensorGet(sensor,'wave');
f = sensorGet(sensor,'spectral qe');
plot(w,f);

%% End
