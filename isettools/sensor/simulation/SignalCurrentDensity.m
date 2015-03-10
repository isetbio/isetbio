function scdImage = SignalCurrentDensity(OI,ISA)
% Estimate the signal current density (current/meter^2) across the sensor surface 
%
%       scdImage = SignalCurrentDensity(OI,ISA)
%
%  This image has a spatial sampling density equal to the spatial sampling
%  of the scene and describes the current per meter.
%
%  We perform the calculation two ways, depending on image size. First, we
%  try to calculate using a quick matrix multiplication. If failed, we try
%  to compute the signal current by looping over all wavebands. It's
%  slower, anyway.
%
% Computational steps:
%
%   The irradiance image in photons (quanta) is multiplied by the spectral
%   QE information. 
%
%   The calculation treats the input data as photons, estimates the
%   fraction of these that are effective, and then turns this into a charge
%   per unit area.  
%
%   Subsequent calculations account for the photodetector area.
%
%   There are many comments in the code explaining each step.
%
% Copyright ImagEval Consultants, LLC, 2003.

q = vcConstants('q');       % Charge per electron

% Hack.  But if we use sceneGet, we get all the data back.
if ~checkfields(OI, 'data', 'photons')   
    warning('Optical image irradiance in photons is required.'); 
    signalCurrentDensityImage = []; %#ok<NASGU>
    return; 
end

% Optical image variables.
oiWaveBinwidth = oiGet(OI,'binwidth');
nRows   = oiGet(OI,'rows');
nCols   = oiGet(OI,'cols');
oiWave  = oiGet(OI,'wave');
oiNWave = oiGet(OI,'nwave');

% Sensor variables
nFilters = sensorGet(ISA,'nfilters');
spectralQE = sensorGet(ISA,'spectralqe'); 
sensorWave = sensorGet(ISA,'wave');

% It is possible that the sensor spectral QE is not specified at the
% same wavelength sampling resolution as the irradiance.  In that case,
% we resample to the lower wavelength sampling resolution.
if ~isequal(oiWave,sensorWave) 
    % Adjust the sensor spectral QE wavelength sampling, in all of the
    % sensor color channels,  to match the irradiance wavelength sampling.
    % We do not change the sensor wavelength data here.
    if length(sensorWave) > 1
        spectralQE = interp1(sensorWave,spectralQE,oiWave,'linear',0);
    elseif ~isequal(sensorWave,oiWave)
        errordlg('Mis-match in sensor and oi wavelength functions.');
    end
end

%  At this point, the spectral quantum efficiency is defined over
%  wavelength bins of size oiWaveBinWidth. To count the number photons in
%  the entire bin, we must multiply by the bin width.
sQE = spectralQE*oiWaveBinwidth;

% Sensor etendue:  In all ISET calculations we treat the etendue (i.e. the
% pixel vignetting) as if it is wavelength independent.  This is an OK
% approximation.  But if we ever want to treat etendue as a function of
% wavelength, we will have to account for it at this point, before we
% collapse all the wavelength information into a single number (the signal
% current density).
%
% If we do that, we may need a space-varying wavelength calculation.  That
% would be computationally expensive.  We aren't yet ready for that level
% of detail.
%
% At present the etendue calculation is incorporated as a single scale
% factor at each pixel and incorporated in the sensorCompute routine. 

% sQE is a wavelength x nSensor matrix, and it includes a conversion
% factor that will maps the electrons per square meter into amps per
% square meter

% Multiply the optical image with the photodetector QE and the color
% filters.  Accumulated this way, we form a current density image at every
% position for all the color filters.
% Output units: [A/m^2]

try
    % This is probably much faster.  But if we are trying to limit the
    % memory size, we should use the other part of the loop that calculates
    % one waveband at a time.
    irradiance = oiGet(OI,'photons');       % quanta/m2/nm/sec
    irradiance = RGB2XWFormat(irradiance);
    
    scdImage =  irradiance * sQE; % (quanta/m2/sec)
    scdImage = XW2RGBFormat(scdImage,nRows,nCols);
    % At this point, if we multiply by the photodetector area and the
    % integration time, that gives us the number of electrons at a pixel.
catch
    % For large images, don't take all of the data out at once.  Do it a
    % waveband at a time.
    scdImage = zeros(nRows,nCols,nFilters);
    for ii=1:oiNWave
        irradiance = oiGet(OI,'photons',oiWave(ii));
        
        for jj=1:nFilters
            scdImage(:,:,jj) = scdImage(:,:,jj) + irradiance*sQE(ii,jj);
        end
    end
end

% Convert the photons into a charge using the constant that defines
% charge/electron.  This is the signal current density (scd) image 
% It has units of quanta/m2/sec/bin * charge/quanta = charge/m2/sec/bin
scdImage = scdImage * q;  

end