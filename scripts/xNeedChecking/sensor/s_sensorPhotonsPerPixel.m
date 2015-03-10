%% s_PhotonsPerPixel
%
% Are you curious about how many total photons are incident on the photodetector
% in a pixel?  Then follow the procedures in this script.
%
% Copyright ImagEval Consultants, LLC, 2005.

% First, use the GUI to make a uniform image (scene window).  

scene  = sceneCreate('uniformee');
oi     = oiCreate;
optics = opticsCreate('human');
oi     = oiSet(oi,'optics',optics);
oi     = oiCompute(scene,oi);

% Then compute
% the irradiance in the optical image window.  Here, we read the irradiance
% in photons.
irradiance = oiGet(oi,'photons');
nWave      = oiGet(oi,'nwave');
wave       = oiGet(oi,'wave');

% Next, set up the ISA and PIXEL data in the sensor window.  We read 
% the size of the photodetector area like this.
sensor = sensorCreate;
pixel = sensorGet(sensor,'pixel');
photonFluxPerPixel = irradiance*pixelGet(pixel,'pdarea');

% Finally, we will look at only a few pixels in the center of the image. 
m = [];
for ii=1:nWave
    m(:,:,ii) = getMiddleMatrix(photonFluxPerPixel(:,:,ii),3);
end

% We calculate the mean spectral power distribution, in photons, by
% converting the RGB image (m) into a matrix that has spatial position
% along the rows and wavelength along the columns (XW format, that is
% space-wavelenth format).
xw = RGB2XWFormat(m);
meanSPD = mean(xw);

% Let's open up a graph window and plot the function
figNum = vcNewGraphWin;
plot(wave,meanSPD);
grid on
set(gca,'ylim',[.95*min(meanSPD),1.05*max(meanSPD)])

% At the top of the window, we print the total number of photons summed
% across all wavelengths as well as the photodetector width.
pdWidth = pixelGet(pixel,'pdWidth','um'); 
txt = sprintf('Total photons: %.0f -- pd width %.1f (um) -- Scene Lum %.0f (cd/m2)',...
    sum(meanSPD),pdWidth,sceneGet(scene,'meanluminance'));
title(txt)
xlabel('Wavelength (nm)');
ylabel('Photons/pixel/nm/sec')

%% End
