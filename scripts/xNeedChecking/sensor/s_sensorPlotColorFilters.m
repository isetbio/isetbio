%% s_plotSensorColorFilters
%
% Plots the transmissivities of color filters in the file stored in the
% ISET data/sensor area. At some point, all of these filter files will be
% moved into data/sensor/colorfilters but some are still in the wrong
% directory.
%
% There is also a routine to create Gaussian color filters with a
% specification of wavelength, wavelength center and bandwidth.  
%
%  cfType = 'gaussian'; wave = [350:850];
%  cPos = 450:50:750; width = ones(size(cPos))*25;
%  fData = sensorColorFilter(cfType,wave, cPos, width);
%  plot(wave,fData)
%
% Copyright ImagEval Consultants, LLC, 2010


%% Look through some files
%
% 'interleavedRGBW' - Nikon RGB and very wide, very sensitive white
% 'gaussianBGRWwithIR' - Gaussian at 450, 580 640 100 nm or so and broad
%        Gaussian at 580 (300 nm)
% 'sixChannel'  - RGBCYM
% 'RGBW' - RGB, limited to 700nm constant W = constant at 0.25 over visible
% 'NikonD1' - Needs filter names
% 'NikonD70' - Needs filter names
% 'NikonD100' - 
% 'NikonD200IR' - Has IR sensitivity 
% 'interleavedRGBW'
% 'M' - magenta
% 'C' - cyan
% 'Y' - yellow
% 'R' - red
% 'G' - green
% 'W' - white
% 'GRBC' - green, red, blue, cyan
%
% 

%
wavelength = 400:1000;

%%
[data,newFilterNames] = ieReadColorFilter(wavelength,'GRBC');
plot(wavelength,data)
xlabel('nm')


