function outSensor = sensorCompute(sensor,oi,showBar)
%Compute sensor response using sensor parameters and optical image data
%
%   inSensor = sensorCompute([sensor],[oi],[showBar = 1])
%
%  This  top-level function  combines the parameters of an image sensor
%  array (sensor) and an optical image (oi) to produce the sensor volts
%  (electrons).
%
% Inputs
%    sensor:   Image sensor, possibly an array of sensors
%    oi:       Optical irradiance data
%    showBar:  Show waitbar if 1, default is off (0)
%
%  On return the voltage data are stored in the sensor.
%
% Return
%    outSensor: Image sensor containing the data.  Possibly an array.
%
%  The computation checks a variety of parameters and flags in the sensor
%  structure to perform the calculation.  These parameters and flags can be
%  set either through the graphical user interface (sensorImageWindow) or
%  by scripts.
%
%  One of the most important flags is the noise flag.  In general, the
%  sensorCompute includes photon and electrical noise in the computation.
%  It can be useful to run the calculation without any noise, or with
%  photon noise only.  To achieve this set the flag as follows:
%
%    sensor = sensorSet(sensor,'noise flag',0);   % No noise
%    sensor = sensorSet(sensor,'noise flag',1);   % Photon only
%    sensor = sensorSet(sensor,'noise flag',0);   % Photon and electrical
%
%
% COMPUTATIONAL OUTLINE:
%
%   This routine provides an overview of the algorithms.  The specific
%   algorithms are described in the routines themselves.
%
%   1. Handle exposure: autoExposure (or not).
%   2. Compute the mean image: sensorComputeImage()
%   3. Etendue calculation
%   4. Noise, analog gain, clipping, quantization
%   5. Correlated double-sampling
%   6. Handle macbeth management
%
%  The value of showBar determines whether the waitbar is displayed to
%  indicate progress during the computation.
%
% See also:  sensorComputeNoise, sensorAddNoise
%
% Key computational flags:
%
%   Exposure Duration
%   Noise calculations
%
% Examples:
%   sensor = sensorCompute;   % Use selected sensor and oi
%   tmp = sensorCompute(vcGetObject('sensor'),vcGetObject('oi'),0);
%
%  Or, compute with specific sensors
%   scene = sceneCreate; scene = sceneSet(scene,'hfov',4);
%   oi = oiCreate; sensor = sensorCreate;
%   oi = oiCompute(oi,scene); sensor = sensorCompute(sensor,oi);
%   vcAddAndSelectObject(sensor); sensorWindow('scale',1);
%
% Copyright ImagEval Consultants, LLC, 2011

%% Define and initialize parameters
if notDefined('sensor'),  sensor = vcGetSelectedObject('sensor'); end
if notDefined('oi'),      oi = vcGetSelectedObject('oi');         end
if notDefined('showBar'), showBar = ieSessionGet('waitbar');      end

wBar = [];

% We allow sensor arrays.  This was necessary as a temporary edit to keep
% the code similar for a while.  Later, I will use thisSensor and simplify
% the logic.
masterSensor = sensor;
clear sensor;

for ss=1:length(masterSensor)   % Number of sensors
    sensor = masterSensor(ss);
    %% Standard compute path
    if showBar, wBar = waitbar(0,sprintf('Sensor %d image:  ',ss)); end
    
    % Determine the exposure model
    integrationTime = sensorGet(sensor, 'integration Time');
    pattern = sensorGet(sensor,'pattern');
    if numel(integrationTime) == 1 && ...
            ( (integrationTime == 0) || sensorGet(sensor,'auto exposure') )
        % The autoexposure will need to work for the cases of 1 value for
        % the whole array and it will need to work for the case in which
        % the exposure times have the same shape as the pattern.  If
        % neither holds then we have to use the vector of numbers that are
        % sent in. We could decide that if autoexposure is on and there is
        % a vector of values we replace them with a single value.
        if showBar, wBar = waitbar(0,wBar,'Sensor: Auto Exposure'); end
        sensor.integrationTime  = autoExposure(oi,sensor);
        
    elseif isvector(integrationTime)
        % We are in bracketing mode, do nothing.
        
    elseif isequal( size(integrationTime),size(pattern) )
        % Find best exposure for each color filter separately
        if sensorGet(sensor,'autoexposure')
            sensor = sensorSet(sensor, 'exp time', ...
                            autoExposure(oi,sensor,[],'cfa'));
        end
    end
    
    %% Calculate current
    % This factor converts pixel current to volts for this integration time
    % The conversion units are
    %
    %   sec * (V/e) * (e/charge) = sec * V / charge = V / current.
    %
    % Given the basic rule V = IR, k is effectively a measure of resistance
    % that converts current into volts given the exposure duration.
    %
    % We handle the case in which the integration time is a vector or
    % matrix, by creating a matching conversion from current to volts.
    q = vcConstants('q');     %Charge/electron
    pixel = sensorGet(sensor,'pixel');
    
    cur2volt = sensorGet(sensor,'integrationTime') * ...
                            pixelGet(pixel,'conversionGain') / q;
    cur2volt = cur2volt(:);
    
    % Calculate the signal current assuming cur2volt = 1;
    if showBar, unitSigCurrent = signalCurrent(oi,sensor,wBar);
    else            unitSigCurrent = signalCurrent(oi,sensor);
    end
    
    
    %% Convert to volts
    % Handle multiple exposure value case.
    if numel(cur2volt) == 1
        % There is only one exposure time.  Conventional calculation
        voltImage = cur2volt*unitSigCurrent;
    else
        % Multiple exposure times, so we copy the unit term into multiple
        % dimensions
        voltImage = repmat(unitSigCurrent,[1 1 length(cur2volt)]);
        % Multiply each dimension by its own scale factor
        for ii=1:length(cur2volt)
            voltImage (:,:,ii) = cur2volt(ii) * voltImage (:,:,ii);
        end
    end
    
    %% Calculate etendue from pixel vignetting
    % We want an array (wavelength-independent) of scale factors that account
    % for the loss of light() at each pixel as a function of the chief ray
    % angle. This method only works for wavelength-independent relative
    % illumination. See notes in signalCurrentDensity for another approach that
    % we might use some day, say in ISET-3.0
    if showBar, waitbar(0.4,wBar,'Sensor image: Optical Efficiency'); end
    sensor    = sensorVignetting(sensor);
    etendue   = sensorGet(sensor,'sensorEtendue');  
    % vcNewGraphWin; imagesc(etendue)
    
    voltImage = bsxfun(@times, voltImage, etendue);
    % vcNewGraphWin; imagesc(voltImage); colormap(gray)
    
    % This can be a single matrix or it can a volume of data. We should
    % consider whether we want to place the volume elsewhere and always
    % have voltImage be a matrix is displayed into the sensor window
    voltImage(voltImage < 0) = 0;
    sensor = sensorSet(sensor, 'volts', voltImage);
    
    % Something went wrong.  Return data empty,  including the noise images.
    if isempty(voltImage),
        % Something went wrong. Clean up the mess and return control to the
        % main processes.
        delete(wBar); return;
    end
    
    %% We have the mean image computed.  We add noise, clip and quantize
    
    % If there is some request for noise (noiseFlag not 0) run the noise.
    if sensorGet(sensor, 'noiseFlag')
        % If you have the mean and all you want is noise, you can just run
        % this.
        sensor = sensorComputeNoise(sensor,wBar);
    end
    
    
    %% Correlated double sampling
    if  sensorGet(sensor,'cds')
        % Read a zero integration time image that we will subtract from the
        % simulated image.  This removes much of the effect of dsnu.
        integrationTime = sensorGet(sensor,'integration time');
        sensor = sensorSet(sensor,'integration time',0);
        
        if showBar, waitbar(0.6,wBar,'Sensor image: CDS'); end
        cdsVolts = sensorComputeImage(oi,sensor);    %THIS WILL BREAK!!!!
        sensor = sensorSet(sensor,'integration time',integrationTime);
        sensor = sensorSet(sensor,'volts', ieClip(sensor.data.volts - cdsVolts,0,[]));
    end
    
    if isempty(sensorGet(sensor,'volts')),
        % Something went wrong.  Clean up the mess and return control to
        % the main processes.
        delete(wBar); return;
    end
    
    %% Macbeth chart management
    % Possible overlay showing center of Macbeth chart
    sensor = sensorSet(sensor,'mccRectHandles',[]);
    
    if showBar, close(wBar); end
    % The sensor structure has new fields at this point, so reassigning to
    % the input sensor array doesn't work.
    outSensor(ss) = sensor;
    
end

end