function sensor = sensorAddNoise(sensor)
% Add electrical and photon noise to the sensor voltage image
%
%    sensor = sensorAddNoise(sensor)
%
% Typically, the sensor structure contains the mean voltage image (without
% noise). Here, we compute the photon noise, sensor electrical noise, and
% quantization error here and add them into the voltage image.
%
% Certain noise terms (fixed pattern noise, seed parameters?) are stored
% and returned in the sensor structure that is returned.  Hence if you want
% to run exactly the same noise simulation again, you can do it by calling
% the function with these seed and related parameters.
%
% An important reason for this routine is to let us generate multiple
% (noisy) samples of the same mean image.
%
% The order of operations are:
%    (a) Add the dark current
%    (b) Add shot noise
%    (c) Add read and reset noise
%    (d) Fixed pattern noise (equivalent to noiseFPN.m)
%    (e) Column fixed pattern noise
%
% See also:  sensorComputeNoise
%
% Copyright ImagEval Consultants, LLC, 2011.

%% Add noise

% We create a noise parameter structure
pixel = sensorGet(sensor,'pixel');

% 0 means do all the noise in the usual way
% 1 means exclude electronic noise, but include shot (photon) noise
% 2 means no electronic noise or shot noise.  Why are you even here?
% 3 could mean do electronic but no shot noise.  But that seems insane.
noiseFlag = sensorGet(sensor,'noise Flag');

% Random noise generator seed issues must be handled here.  This is tough
% to make sure we can absolutely replicate the noise to test code accuracy.
if sensorGet(sensor, 'reuse noise')
    % Get the state stored in the sensor
    noiseSeed = sensorGet(sensor,'noise seed');
    if isempty(noiseSeed)
        noiseSeed = rng;
        sensorSet(sensor, 'noise seed', noiseSeed);
    else
        % Initialize with current seed.
        rng(noiseSeed);
    end
    
    % We should really get rid of the stored dsnu/prnu images. Although
    % we might decide we want it transiently so we can plot it or check it.
    % So it is still in the structure, but not used in any computational
    % path.
    sensor = sensorSet(sensor, 'dsnu image', []);
    sensor = sensorSet(sensor, 'prnu image', []);
else
    % Not reusing.  But remember the initial noise state for this
    % calculation
    noiseSeed = rng;
    sensor = sensorSet(sensor, 'noise seed', noiseSeed);
end

%% Perform the noise addition steps here
nExposures = sensorGet(sensor, 'nExposures');
eTimes = sensorGet(sensor,'exposure times');
volts = sensorGet(sensor,'volts');

% Check if nExposures and the data are consistant
if nExposures ~= size(volts, 3)
    disp('For multisamples in human sensor, use coneAbsorptions instead');
    error('nExposures and data mismatch');
end

for ii=1:nExposures 
    vImage = volts(:,:,ii);
    
    % Add the dark current At this point the noise dark current is the same
    % at all pixels. Later, we apply the PRNU gain factor to the sum of the
    % signal and noise, so that the noise dark current effectively varies
    % across pixels.  Sam Kavusi says that this variation in gain (also
    % called PRNU) is not precisely the same for signal and noise.  But we
    % have no way to assess this for most cases, so we treat the PRNU for
    % noise and signal as the same until forced to do it otherwise.
    if noiseFlag > 1
        vImage = vImage + pixelGet(pixel,'dark Voltage') * eTimes(ii);
        sensor = sensorSet(sensor,'volts',vImage);
    end
    
    % Add shot noise  (equivalent to noiseShot.m)
    % Always do this.  It must be done after adding the dark current.
    % For an ideal sensor, there is no dark current so this calculation is
    % simply photon noise.
    if  noiseFlag > 0
        vImage = noiseShot(sensor);
    end
    
    % Add read noise  (equivalent to noiseRead.m)
    if noiseFlag > 1
        readVolts = pixelGet(pixel,'read noise volts');
        vImage = vImage + readVolts * randn(size(vImage));
        sensor = sensorSet(sensor,'volts',vImage);
        
        % Fixed pattern noise (equivalent to noiseFPN.m)
        vImage = noiseFPN(sensor);
        sensor = sensorSet(sensor,'volts',vImage);
        
        % Column fixed pattern noise (equivalent to noiseColumnFPN.m)
        vImage = noiseColumnFPN(sensor);
    end
    
    % That's is.  Store it in the volume image ...
    volts(:,:,ii) = vImage;
end

sensor = sensorSet(sensor,'volts',volts);

end