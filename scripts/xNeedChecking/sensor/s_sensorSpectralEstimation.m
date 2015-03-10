%% s_sensorSpectral
%
%  Illustrate an experimental approach to characterizing sensor color
%  responsivity
%   
% Copyright ImagEval Consultants, LLC, 2010.

%% ISET, ieMainW('visible','off')

%% Make a uniform scene 
scene = sceneCreate('uniformee');
wave  = sceneGet(scene,'wave');

oi     = oiCreate('default',[],[],0);
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'optics model','skip');
oi     = oiSet(oi,'optics',optics);

sensor = sensorCreate; 
sensor = sensorSet(sensor,'size',[64 64]);
sensor = sensorSet(sensor,'auto exposure',1);

% You can experiment using color filters with Gaussian transmissivity, too
% cfType = 'gaussian';
% cPos = [450,530,600]; width = [40,40,40];
% fData = sensorColorFilter(cfType,wave, cPos, width);
% fData = fData(:,[3 2 1]);
% sensor = sensorSet(sensor,'color filters',fData);
% plot(wave,fData(:,1),'-r',wave,fData(:,2),'-g',wave,fData(:,3),'-b')
% sensorGet(sensor,'filter names')

% The scene must always be larger than the sensor field of view.
scene = sceneSet(scene,'fov',sensorGet(sensor,'fov',scene,oi)*1.5);

%% Generate SPDs to use in test scene
waveStep = 50;
cPos     = (wave(1):waveStep:wave(end));  % Center wavelengths
widths   = waveStep/2;
nLights = length(cPos);

% These are Gaussian shaped SPDs.  cPos is the center position and widths
% is the width of the Gaussian.  You can plot them below.
spd = zeros(length(wave),nLights);
for ii=1:nLights
    spd(:,ii) = exp(-1/2*( (wave-cPos(ii))/(widths)).^2);
end
spd = spd*10^16;  % Make them a reasonable number

% vcNewGraphWin; plot(wave,spd)

%% Create one of the spectral scenes, compute the oi and the sensor, save
% the data
eTime        = zeros(1,nLights);
nFilters     = sensorGet(sensor,'nfilters');
volts        = cell(1,nFilters);
responsivity = zeros(nFilters,nLights);

for ii=1:nLights
    
    % Make a scene with a particular spectral power distribution (spd). The
    % code has to arrange the data into the proper 3d matrix format.
    spdImage = repmat(spd(:,ii),[1,32,32]); 
    spdImage = permute(spdImage,[2 3 1]);
    scene    = sceneSet(scene,'cphotons',spdImage);
    % vcAddAndSelectObject(scene); sceneWindow;
    
    % Compute the optical image
    oi     = oiCompute(scene,oi);
    vcAddAndSelectObject(oi); % oiWindow;
    
    % Compute the sensor response.
    sensor    = sensorCompute(sensor,oi,0);
    eTime(ii) = sensorGet(sensor,'Exposure Time','sec');
    % vcAddAndSelectObject(sensor); sensorImageWindow;

    % Calculate volts/sec for each of the channels at this wavelength
    for jj=1:nFilters
        volts{jj} = sensorGet(sensor,'volts',jj);
        responsivity(jj,ii) = mean(volts{jj})/eTime(ii);  % volts/sec
    end
end

% plot(wave,sensorGet(sensor,'color filters'))

%% Estimate filters from the measurements
%
% Use  linear estimation to calculate filters from responsivities
%
%    responsivity = filters*spd;
% 
% So figure that filters are weighted sums of the spd's
%
%   filters = wgt*spd'
%
% Then, wgt = responsivity*inv(spd'*spd);
%       filters = wgt*spd';
%
wgt = responsivity/(spd'*spd);  % Solve for weights
cFilters = (wgt*spd')';         % Solve for filters

% Normalize to peak of 1
cFilters = cFilters/max(cFilters(:));

% The estimates should match the sensor color filters
f = sensorGet(sensor,'color filters');
f = f/max(f(:));
vcNewGraphWin;
subplot(1,2,1), plot(wave,cFilters); grid on; set(gca,'ylim',[0 1])
title('Estimate'); xlabel('Wavelength (nm)'); ylabel('Spectral QE')
subplot(1,2,2), plot(wave,f); grid on; set(gca,'ylim',[0 1])
title('Sensor');


%% End
