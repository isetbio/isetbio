% v_wvfPlot
%
%  Illustrate ways to create plots of the wvf structure using the wvfPlot
%  call.
%
% (BW) Wavefront Toolbox Team, 2014
%
%%
s_initISET

%%  Start with a clean structure
wvf = wvfCreate;
wave = 550; wvf = wvfSet(wvf,'calc wave',wave);
wvf = wvfComputePSF(wvf);

% Make the plot in microns
unit = 'um';
[u,p]= wvfPlot(wvf,'1d psf space',unit,wave);
set(p,'color','k','linewidth',2)
u

% Normalize the plot
unit = 'um';
[u,p]= wvfPlot(wvf,'1d psf space normalized',unit,wave);
set(p,'color','b','linewidth',2)
u

% Make the plot distance axis millimeters
unit = 'mm';
[u,p]= wvfPlot(wvf,'1d psf space normalized',unit,wave);
set(p,'color','r','linewidth',3,'linestyle',':')

%% Show the return arguments

% data values
u

% Figure properties that can be set.
get(p)


%%  Change the calculated PSF wavelength and plot again
wave = 460; wvf = wvfSet(wvf,'calc wave',wave);
wvf = wvfComputePSF(wvf);
unit = 'min';
wvfPlot(wvf,'image psf angle',unit,wave);

%% A multiple axis window

vcNewGraphWin([],'tall');
subplot(3,1,1), wvfPlot(wvf,'1d psf space',unit,wave,'no window');
subplot(3,1,2), wvfPlot(wvf,'1d psf space normalized',unit,wave,'no window');
subplot(3,1,3), wvfPlot(wvf,'image psf','um',wave,20,'no window');

%% Pupil phase
unit = 'mm'; maxMM = 2;
wvfPlot(wvf,'image pupil phase',unit,wave,maxMM);

wvfPlot(wvf,'image pupil amp',unit,wave,maxMM);

%%  Mesh plots of the psf in angle and space

unit = 'min'; maxMIN = 10;
wvfPlot(wvf,'2d psf angle',unit,wave,maxMIN);

unit = 'mm'; maxMM = .050;
wvfPlot(wvf,'2d psf space',unit,wave,maxMM);

% These are linepairs / unit and maximum frequency
unit = 'mm'; maxF = 300;
wvfPlot(wvf,'2d otf',unit,wave,maxF);

%% END