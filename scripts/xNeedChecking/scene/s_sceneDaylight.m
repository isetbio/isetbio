%% s_sceneDaylight
%
%  Create a daylight illuminant with a specific correlated color
%  temperature, or just by choosing random weights for the daylight bases
%  in the CIE standard daylight representation.
%
% Copyright ImagEval Consultants, LLC, 2010.


%% Create a few daylight spectra

wave = 400:770;                % Wavelength in nanometers
cct = 4000:1000:10000;  % Correlated color temperature
spd = daylight(wave,cct,'photons');

% Calculate the luminance of these SPDs
lum = ieLuminanceFromPhotons(spd',wave(:));

% Scale the luminance to 100 cd/m^2
spd = spd*diag(100./lum);

vcNewGraphWin
plot(wave,spd);
grid on
xlabel('Wavelength')
ylabel('Photons (q/sr/m^2/s)');

%% Same data in energy units
spd = daylight(wave,cct,'energy');

% Calculate the luminance of these SPDs
lum = ieLuminanceFromEnergy(spd',wave(:));

% Scale them to 100 cd/m2
spd = spd*diag(100./lum);

vcNewGraphWin
plot(wave,spd);
grid on
xlabel('Wavelength')
ylabel('Energy (watts/sr/m^2/s)');

%% Plot the CIE daylight basis functions
dayBasis = ieReadSpectra('cieDaylightBasis',wave);

vcNewGraphWin
p = plot(wave,dayBasis,'k-');
set(p,'linewidth',2);
xlabel('Wavelength')
ylabel('Energy (relative)');
grid on

%% An example of three daylights 

% The mean, and +/1 the first coefficient.
wgts = [1 0 0 ; 1, 1, 0; 1 -1 0]';
vcNewGraphWin
plot(wave,dayBasis*wgts)
grid on;
xlabel('Wavelength')
ylabel('Energy (relative)');

