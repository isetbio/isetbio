%% s_sensorSpatialNoiseDSNU
%
%   Illustrate an experimental approach to measuring DSNU .
%
% The method is described in the paper by Farrell, Feng and Kavusi. We
% simulate a large number of short exposure durations to a black scene. We
% average across the multiple images to eliminate the (presumably unbiased)
% read noise.  We estimate of the DSNU as the standard deviation across the
% pixels.
%   
% Copyright ImagEval Consultants, LLC, 2010.


%% Make a black scene and a sensor
scene = sceneCreate('uniformee');
darkScene = sceneAdjustLuminance(scene,1e-8);

oi = oiCreate('default',[],[],0);

sensor = sensorCreate; 
sensor = sensorSet(sensor,'size',[196 196]);

% The scene must always be larger than the sensor field of view.
darkScene = sceneSet(darkScene,'fov',sensorGet(sensor,'fov')*1.5);

%% Compute the optical image
darkOI = oiCompute(darkScene,oi);

%% Set sensor parameters
dsnuLevel = 0.05;       % Std. dev. of offset in volts
prnuLevel = 0.2;      % Std. dev of gain, around 1, as a percentage
readNoise = 0.001;      % Read noise in volts

% Set a brief exposure time for DSNU estimation.  Because of the short
% exposure, the PRNU level is irrelevant.  The read noise can matter.  If
% it is very large, you must average over more trials.
expTime = 0.001; 

sensor = sensorSet(sensor,'DSNU level',dsnuLevel);
sensor = sensorSet(sensor,'PRNU level',prnuLevel);
pixel  = sensorGet(sensor,'pixel');
pixel  = pixelSet(pixel,'Read noise volts',readNoise);
sensor = sensorSet(sensor,'pixel',pixel);
sensor = sensorSet(sensor,'Exposure Time',expTime);

% How many color filters?  Normally 3 and we use the 2nd one in a Bayer.
% But sometimes we might run this script with a monochrome.
nFilters = sensorGet(sensor,'nfilters');

%%  Acquire multiple short exposures of the dark image
clear volts

% We take the image multiple times so we can average out the read noise
nRepeats = 25;

wBar = waitbar(0,'Acquiring images');
nSamp = prod(sensorGet(sensor,'size'))/2;
volts = zeros(nSamp,nRepeats);

for ii=1:nRepeats
    waitbar(ii/nRepeats,wBar);
    sensor = sensorCompute(sensor,darkOI,0);
    if nFilters == 3
        volts(:,ii) = sensorGet(sensor,'volts',2);
    elseif nFilters == 1
        tmp = sensorGet(sensor,'volts');
        volts(:,ii) = tmp(:);
    end
end
close(wBar);

% vcAddAndSelectObject(sensor); sensorImageWindow;

%% Estimate the standard deviation across the sensor

% If we average across the repeated measures, we get a number that
% describes the mean offset.
meanOffset = mean(volts,2);

% The histogram gives a sense of the distribution of these offsets, which
% is the DSNU
vcNewGraphWin; 
hist(meanOffset,50)
grid on;
title('Offsets averaged across reads')

% Notice that the values are truncated at 0 because the voltages are
% clipped at 0.  So, to estimate the true DSNU we multiply the std of the
% offsets by 2.  This is not a very good idea.  Rather, we should take the
% distribution of values above 0 and fit it with a truncated Gaussian.
% That Gaussian should have close to a zero mean and we can use the
% standard deviation as a DSNU.

% But I don't want to write an fminsearch routine for this now.  So, as a
% brief hack we just flip the positive values negative and recompute the
% standard deviation.  It's not right, but it's easy for this tutorial.

% Find the values that are away from 0, and assume that we would give a
% symmetric set of negative values.  Compute the std of that joint set.
z = (meanOffset ~= 0.0);
m1 = meanOffset(z); m1 = m1(:)';
estimatedDSNU = std([-m1,m1]);

fprintf('DSNU:  Estimated: %.3f and Set %.3f\n',estimatedDSNU,dsnuLevel);

%% End
