%% s_sensorStackedPixels
%
% Illustrates how to simulate a sensor with stacked pixels, as in the TFD
% and Foveon sensors.  Comments and descriptions are inserted below.
%
% Copyright ImagEval Consultants, LLC, 2012.

%% Initialize ISET
s_initISET

%% Initialize a simple scene

%Following value specifies horizontal field of view code in degrees.  A
%good value is approximately linearly related to the size of the images.
%If the field of view is too large, the output image will contain more
%pixels than the input scene.  If the field of view is too small, the
%output image will be much smaller than the original scene.  A
%precautionary error will occur if the field of view is too large or small.
horizontalFOV = 4;

% Illuminant filenames (should match a file on the path and contain wavelength and data)
% inputilluminantname ='D65';      %illuminant that the scene will be under
% outputilluminantname='D65';     %illuminant that you want the output to be rendered under

% Ideal color bands to be estimated name of file on path that defines
% desired output bands (contains wavelength and data).  There could be a
% function for defining the spectral curves for different assumptions about
% the FOVEON or TFD designs.  These different filter files would be read in
% here.
%
% For examples on how to create such filter files see xyzQuantaCreate.m, or
% have a look at ieSaveSpectralFile.m
filterFile = 'XYZQuanta';

%Following is to specify the scene mean luminance 
meanLuminance = 100;

%% Get the scene - see script un Data/Scenes
scene = sceneCreate;

% Set the scene mean luminance somewhere reasonable
scene = sceneAdjustLuminance(scene,meanLuminance);
scene = sceneSet(scene,'hfov',horizontalFOV);
wave = sceneGet(scene,'wave');

% vcAddAndSelectObject(scene); sceneWindow;

%% Build the OI
oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'f number',4);
optics = opticsSet(optics,'focal length',3e-3);   % units are meters

%% Compute optical image
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); %oiWindow;
%% Simulated sensor

% Create a monochrome sensor.  We will reuse this structure to compute each
% of the complete color filters.
sensor = sensorCreate('monochrome');
sensor = sensorSet(sensor,'wave',wave);

% Set sensor parameters
sensor = sensorSet(sensor,'filter spectra',ones(length(wave),1));
sensor = sensorSet(sensor,'irfilter',ones(length(wave),1));
sensor = sensorSet(sensor,'quantizationMethod','analog');

% Set pixel parameters - could pull these parameters out as variables.
% These are just examples.
pixel = sensorGet(sensor,'pixel');
pixel = pixelSet(pixel,'spectralQE',ones(length(wave),1));
pixel = pixelSet(pixel,'size',[2.2e-6 2.2e-6]);                % Pixel Size in meters
pixel = pixelSet(pixel,'conversion gain', 2.0000e-004);        % Volts/e-
pixel = pixelSet(pixel,'voltage swing', 1.8);                  % Volts/e-
pixel = pixelSet(pixel,'dark voltage', 1e-005);                % units are volts/sec
pixel = pixelSet(pixel,'read noise volts', 1.34e-003);         % units are volts

sensor = pixelCenterFillPD(sensor, 0.45);

sensor = sensorSet(sensor,'pixel',pixel);
sensor = sensorSet(sensor,'dsnu level',14.1e-004); % units are volts
sensor = sensorSet(sensor,'prnu level',0.002218);  % units are percent

%Following sets horizontal field of view to desired value
% Maybe we should put in a flag to allow for vertical and horizontal.
% sceneGet(scene,'vfov')
% sceneGet(scene,'hfov')
sensor = sensorSetSizeToFOV(sensor,horizontalFOV,scene,oi);   % deg of visual angle
sensorsize = sensorGet(sensor,'size');
rows = sensorsize(1);
cols = sensorsize(2);
if isodd(rows), rows = rows+1; end
if isodd(cols), cols = cols+1; end
sensor = sensorSet(sensor,'size',[rows cols]);

%% Read the CIE XYZ filters for quantal calculations

wave = sensorGet(sensor,'wave');
fSpectra = ieReadSpectra(filterFile,wave);   %load and interpolate filters
filterNames = {'rX','gY','bZ'};
nChannels = size(fSpectra,2); 
% vcNewGraphWin; plot(wave,fSpectra)

%% Loop on the number of filters and calculate full sensor plane values

% For each of the filter transmissivities, compute the monoSensor photons.
sz = sensorGet(sensor,'size');

% We will store the image values here
im = zeros(sz(1),sz(2),nChannels);

for kk=1:nChannels
    s = sensorSet(sensor,'filterspectra',fSpectra(:,kk));
    s = sensorSet(s,'Name',sprintf('Channel-%.0f',kk));
    s = sensorCompute(s,oi,0);
    im(:,:,kk) = sensorGet(s,'volts');
end

%% The multiple channel data (in volts) are stored in the variable im().

% The data in the variable 'im' are the voltages that one would obtain from
% the sensor with stacked pixels.  The experiments that could be run
% involve adjusting the filters and other optics and sensor properties. 
%
% You can visualize the data summed across the channels using:
% hcimage(im), or hcimage(im,'montage')
%
% Below shows howto create a mosaicked version of the data. It is possible
% to create a virtual camera image (vci) from these data without applying
% any demosaicking.  Ask me if you want help with this.  Or perhaps you can
% just see the steps in the function vcimageCompute.

%% If you would would like to create a standard mosaicked sensor 
%  with these data, say with a Bayer CFA pattern, you can do this

% Suppose you want a Bayer mosaic
cfaPattern = [1 2 ; 2 3];

% Create a copy of the sensor and fill it with the filters and cfa pattern
sensor2 = sensor;
sensor2 = sensorSet(sensor2,'cfa pattern',cfaPattern);
sensor2 = sensorSet(sensor2,'filter spectra',fSpectra);
sensor2 = sensorSet(sensor2,'filter names',filterNames);

% This function converts the full image data into the sampled mosaic
sensorMosaic = sensorRGB2Plane(im, cfaPattern);

% Put the data into the new sensor
sensor2 = sensorSet(sensor2,'volts',sensorMosaic);

% Visualize it.
vcAddAndSelectObject(sensor2); sensorImageWindow;

%% To render the data in the image processor window
% vci = vcimageCreate;
% vci = vcimageCompute(vci,sensor2);
% vcAddAndSelectObject(vci); vcimageWindow;

%% End
