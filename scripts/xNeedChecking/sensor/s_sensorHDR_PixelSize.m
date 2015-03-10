% s_sensorHDR_PixelSize
%
% Integrate data from multiple pixel sizes.
%
% The idea of the simulation is to capture the image data with a series of
% sensors, each has a different pixel size. But all the simulated sensors
% have the same dye size, so it is possible to figure out which row/col
% values correspond in the different sensors. 
% 
% An HDR algorithm might extract data data from corresponding positions in
% the different simulated sensors and integrate these values into a single
% image. 
%
% Copyright ImagEval Consultants, LLC, 2010

% It is convenient to start up ISET and hide the main window.  This makes
% it easier to sometimes visualize the results of different computations
ISET
ieMainW('visible','off')

%% Read a high dynamic range scene

fName = fullfile(isetRootPath,'data/images/multispectral/Feng_Office-hdrs.mat');

% Read in a multispectral file with high dynamic range and a mean level of
% 200 cd/m2
[scene,fname] = sceneFromFile(fName,'multispectral',200);

% vcAddAndSelectObject(scene); sceneWindow
%%  Create an optical image with the default lens parameters (f# = 2.0)
oi = oiCreate;
oi = oiCompute(scene,oi);
% vcAddAndSelectObject(oi); oiWindow

%%  A series of sensors with different pixel sizes but same dye size

clear psSize;
pSize = [1 2 4];        % Microns
dyeSizeMicrons = 512;   % Microns
fillFactor = 0.5;       % Percentage of pixel containing photodector

% We will have a cell array of sensors
sensor = cell(length(pSize),1);

% The sensor will be the same as the base sensor, but we will adjust their
% pixel sizes (keeping dye size constant)
baseSensor = sensorCreate('monochrome');             % Initialize
baseSensor = sensorSet(baseSensor,'expTime',0.003);  % 3 ms exposure time

% This is the base processor image. We store the rendered image here 
baseProcessor = vcimageCreate;    
vci = cell(length(pSize),1);

%% Run the main  loop

% We simulate a series of monochrome sensors with different pixel sizes We
% then render the images and place them in the virtual camera image window.
% We can leaf through them and see the effects of scaling the pixel sizes.
for ii=1:length(pSize)
    
    % Adjust the pixel size (meters)
    pixel = sensorGet(baseSensor,'pixel');
    pixel = pixelSet(pixel,'size',[pSize(ii) pSize(ii)]*1e-6);
    
    % Stuff the pixel back into the sensor structure
    sensor{ii} = sensorSet(baseSensor,'pixel',pixel);
    
    % Make sure the photodetector has the proper fill factor
    sensor{ii} = pixelCenterFillPD(sensor{ii},fillFactor);

    %Adjust the sensor row and column size so that the sensor has a constant
    %field of view.
    sensor{ii} = sensorSet(sensor{ii},'rows',round(dyeSizeMicrons/pSize(ii)));
    sensor{ii} = sensorSet(sensor{ii},'cols',round(dyeSizeMicrons/pSize(ii)));

    sensor{ii} = sensorCompute(sensor{ii},oi);
    sensor{ii} = sensorSet(sensor{ii},'name',sprintf('pSize %.1f',pSize(ii)));
    vcAddAndSelectObject(sensor{ii}); 
    % sensorImageWindow;
    
    vci{ii} = vcimageCompute(baseProcessor,sensor{ii});
    vci{ii} = imageSet(vci{ii},'name',sprintf('pSize %.1f',pSize(ii)));
    vcAddAndSelectObject(vci{ii}); 
end


%% Notes - For computational modeling based on the simulated sensor

% Have a look at the series of images that were created. They have
% different spatial resolutions and row/col sizes.
%
% This brings up the window.  Click around. We suggest setting the display
% gamma (text box, lower left of the window) to about 0.6 
vcimageWindow;  

% You can extract the raw data from different sources. To read out the
% voltages from the sensor directly you can use
ii = 2;
v = sensorGet(sensor{ii},'volts');
size(v)
% The size of v is equal to the number of pixels on the sensor

% To read out the rgb data from the virtual camera image (vci) structures
% you can use
ii = 1;
dv = imageGet(vci{ii},'result');
size(dv)
% The dv are RGB images, but since they are acquired with a monochrome
% camera the RGB values are equal.



