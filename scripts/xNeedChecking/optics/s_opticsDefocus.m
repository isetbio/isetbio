% s_opticsDefocus
%
% Illustrate the OTF and linespreads for defocus under different
% conditions.  
%
% Copyright ImagEval Consultants, LLC, 2011

%% Example: Human lens without chromatic aberration or williams factor.
%
sampleSF = 0:60;  % Spatial frequency
optics = opticsCreate('human');
optics = initDefaultSpectrum(optics);
wave = opticsGet(optics,'wave');

% No defocus per wavelength, so not really human yet.
D = zeros(size(wave));             % No wavelength dependent defocus
otf = opticsDefocusCore(optics,sampleSF,D);

vcNewGraphWin; 
mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%% Set defocus to human
wave = opticsGet(optics,'wave');
D = humanWaveDefocus(wave);
otf = opticsDefocusCore(optics,sampleSF,D);

vcNewGraphWin; 
mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%% Bigger pupil, less diffraction
optics = opticsSet(optics,'focal length',2*opticsGet(optics,'focal length'));
opticsGet(optics,'aperture diameter')
otf = opticsDefocusCore(optics,sampleSF,D);
vcNewGraphWin; mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%% No wavelength dependent defocus (diffraction limited)

D = zeros(size(wave));             % No wavelength dependent defocus
otf = opticsDefocusCore(optics,sampleSF,D);
vcNewGraphWin; mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%%  This is the default lens case

optics = opticsCreate;
optics = initDefaultSpectrum(optics);

sampleSF = 0:60;  % Spatial frequency
wave = opticsGet(optics,'wave');

% Defocus
D = zeros(size(wave));             % No wavelength dependent defocus

% Must deal with complex values somehow ...
otf = opticsDefocusCore(optics,sampleSF,D);

vcNewGraphWin; mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%%  Introduce some wavelength-dependent defocus.

% Defocus
D = 0.5*ones(size(wave));             % No wavelength dependent defocus

% We have a problem with the sampleSF = 0 case.  There is some scaling off
% in the defocus calculation.  Here we just scale it up.
% We should probably put this in the opticsDefocusCore routine itself, I am
% afraid.
otf = opticsDefocusCore(optics,sampleSF,D);
% s = otf(1,1); otf = otf/s;

vcNewGraphWin; mesh(sampleSF,wave,otf)
set(gca,'zlim',[0 1]);
xlabel('Spatial Frequency (cyc/deg)'), ylabel('wave'), zlabel('OTF')

%% End