% s_RawMCC2Sensor
%
% This script converts the a tiff file into ISET sensor format.  The TIFF
% file represents sensor data from a GBRG sensor.
%
% This script is intended to illustrate how you can use ISET to start with
% a raw image data from a camera, insert those data into the ISET pipeline,
% and then continue to explore algorithms and processing.
%
% The script loads an MCC simulation (but with slightly off colors).  It
% then uses sensorCCM to find the sensor color conversion matrix.
%
% In this example, the color rendering of the original image is not very
% good.  So the delta E values are large.
%
% Copyright ImagEval Consultants, LLC, 2010.
%% Read the TIF file

% This is an approximation to the Gretag in some sensor
fName = 'gbrgMCCSensor.tif';
mosaic = imread(fName);   % We treat the data in this file as sensor volts

% imtool(mosaic)

%%  For a normal image, we would probably crop the data

% This is the code you could use.  BUt we don't run it here because for
% this demonstration image the sizes and so forth are OK.

% figure(1), clf

% The image is big.  Show an image and take the part you want to analyze.
% In this imcrop routine you drag to select the rectangle.  When you are
% done, right click and select 'crop image'.
% [croppedMosaic,rect] = imcrop(mosaic);
% rect = round(rect);

% Adjust the rect so we fall neatly on the bayer sampling grid.
% We want the (xmin,ymin) values to both be odd. 
% if ~isodd(rect(1)), rect(1)=rect(1)+1; end
% if ~isodd(rect(2)), rect(2)=rect(2)+1; end

% We want the (width,height) values to both be even.  Matlab's imcrop
% basically adds one more pixel than you want.  So, annoyingly, we must
% make the width and height odd, so we get an even number of pixels out. 
% if ~isodd(rect(3)),  rect(3)=rect(3)+1; end
% if ~isodd(rect(4)),  rect(4)=rect(4)+1; end

% Recrop to make sure we have an even number of samples, matching the 2x2
% of a Bayer array.
% croppedMosaic = imcrop(mosaic,rect);

% Have a look
% imtool(croppedMosaic)


%% Create the sensor structure.

% Simulated camera sensor is gbrg
sensor = sensorCreate('bayer (gbrg)');
sensor = sensorSet(sensor,'Name','Sensor demo');

%% Set the sensor properties - not needed for this

% But if you know them you could set them and simulate other things

% sensor = sensorSet(sensor,'integrationTime',0.1);   %100 ms 
% sensor = sensorSet(sensor,'autoexposure',1);
% sensor = sensorSet(sensor,'exposuretime',0.02); % in units of seconds
% sensor = sensorSet(sensor,'rows',960);
% sensor = sensorSet(sensor,'cols',1280);
% sensor = sensorSet(sensor,'dsnulevel',0);  %mV 
% sensor = sensorSet(sensor,'prnulevel',0.007); % specified in %

%% Set the pixel properties - not needed, for this

% But if you know them you could set them and simulate other things

% pixelSize      = 3.75*1e-6;     %Pixel size (m)
% voltageswing   = 1.8;           % Volts
% darkvoltage    = 0 ;            %V/sec/pixel
% readnoisevolts = 0.00131;       % standard deviation in V
% wellcapacity   = 15000;
% conversiongain = voltageswing/wellcapacity ; % Volts/e-
% %end pixel data
% 
% pixel = sensorGet(sensor,'pixel');
% pixel = pixelSet(pixel,'size',[pixelSize pixelSize]);   
% pixel = pixelSet(pixel,'conversiongain',conversiongain);        
% pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;             
% pixel = pixelSet(pixel,'readnoisevolts',readnoisevolts);   
% 
% % Stuff the pixel back into the sensor structure
% sensor = sensorSet(sensor,'pixel',pixel);

%% Attach the volts to the sensor
pixel  = sensorGet(sensor,'pixel');
vSwing = pixelGet(pixel,'voltageSwing',1.8);

% We want to scale the digital values in the croppedMOsaic so that the max
% is equal to the voltage swing and the min is equal to the dark level.  
mn = double(min(mosaic(:)));
mx = double(max(mosaic(:)));
volts = ((double(mosaic) - mn)/(mx - mn))*vSwing; 
% figure; hist(volts(:),50)

sensor = sensorSet(sensor,'size',size(volts));
sensor = sensorSet(sensor,'volts',volts);

%To view image in GUI
vcAddAndSelectObject(sensor);
sensorImageWindow;

%% Interactively determine Color Conversion Matrix CCM

% Click on the outer corners of the patches in the order described in the
% message within the Sensor Window.  When you are done, right click.  This
% is a pretty big image so it will take 5-10 secs to compute.
sensorCCM;

% When you run this line, the matrix from sensorCCM is printed into the
% command window and some metrics are run and graphs are produced.  We
% copied the matrix values below and use them for rendering.

%%  Render the image without the CCM

% First, compute with the default properties.  This uses bilinear
% demosaicing, no color conversion or balancing.  The sensor RGB values are
% simply set to the display RGB values.

% Create a display image with basic attributes
vci = vcimageCreate;
vci = imageSet(vci,'name','No Correction');
vci = imageSet(vci,'scaledisplay',1);
vci = imageSet(vci,'renderGamma',0.6);
vci = vcimageCompute(vci,sensor);
vcAddAndSelectObject(vci);
vcimageWindow;

%% Render the image with the color conversion matrix

% For this data set, the white patch is saturated so the matrix isn't quite
% right.

vci = vcimageCreate;
vci = imageSet(vci,'name','CCM Correction');
vci = imageSet(vci,'scaledisplay',1);
vci = imageSet(vci,'renderGamma',0.6);

% In the sensor window we used the pulldown under 
% Analyze | Color | Color Conversion Matrix 
% to find an optimal transform matrix for the MCC
m = [ ...
   0.9205   -0.1402   -0.1289
   -0.0148    0.8763   -0.0132
   -0.2516   -0.1567    0.6987];
vci = imageSet(vci,'sensor conversion transform',m);

% We set  the other transforms to the identity, so that the product of all
% the transforms is just the one above.
vci = imageSet(vci,'illuminant correction transform',eye(3,3));
vci = imageSet(vci,'ics2Display Transform',eye(3,3));

% We set the vci to not ask any questions, just use the current matrices.
vci = imageSet(vci,'sensor conversion method','current matrix')

% Compute and show.
vci = vcimageCompute(vci,sensor);
vcAddAndSelectObject(vci); vcimageWindow
    
%% End
