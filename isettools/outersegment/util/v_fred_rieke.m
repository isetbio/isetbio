%% v_fred_rieke
%
%    Compare cone isomerization results between ISETBIO and Rieke's paper
%
%  Reference:
%    Charles A Hass, et. al, Chromatic detection from cone photoreceptors
%    to V1 neurons to behavior in rhesus monkeys
%
%  (HJ) ISETBIO TEAM, 2014

%% Create scene
%  create a half gray scene on sony OLED display
wave = (400:10:700)';
fov  = 1;
d = displayCreate('OLED-Sony', wave);
I = 0.75 * ones(100); % 0.75 turns to 0.5 in luminance after gamma table
scene = sceneFromFile(I, 'rgb', [], d);
scene = sceneSet(scene, 'h fov', fov);

%% Create human optics
%  According to the paper, pupil size should be of area 12.6 mm^2 (4 mm
%  pupil diameter)
pupil_size = 4; % 4 mm diameter
oi = oiCreate('wvf human', pupil_size);
oi = oiCompute(scene, oi);

%% Create sensor
% According to the paper, cone collecting area is 0.6 um^2
%  macular pigment transmittancewas scaled to 0.35 at 460 nm
%  lens transmittancewas scaled to 1 at 400 nm
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'wave', wave);

pixel = pixelCreate('human', wave);
pixel = pixelSet(pixel, 'pd width', 0.774e-6); % photo-detector width
pixel = pixelSet(pixel, 'pd height', 0.774e-6);
sensor = sensorSet(sensor, 'pixel', pixel);

% lens transmittance was scaled to 1 at 400 nm
% THIS PART IS NOT RIGHT AT THIS POINT...NOW IT'S TO EQUIVALENT TO IGNORE
% THE LENS EFFECT...
lens = sensorGet(sensor, 'human lens');
lens_trans = lensGet(lens, 'transmittance');
unit_density = lensGet(lens, 'unit density');
scale_factor = log10(lens_trans(wave==400))/unit_density(wave==400);
lens_density = lensGet(lens, 'density');
lens = lensSet(lens, 'density', lens_density + scale_factor);
sensor = sensorSet(sensor, 'human lens', lens);

% macular pigment absorbance was scaled to 0.35 at 460 nm
macular = sensorGet(sensor, 'human macular');
macular = macularSet(macular, 'density', 0.35);
macular = macularSet(macular, 'density', 0);
sensor = sensorSet(sensor, 'human macular', macular);

%% Compute cone absorptions
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', 1);
sensor = sensorComputeNoiseFree(sensor, oi);

photons = gather(sensorGet(sensor, 'photons'));
coneType = sensorGet(sensor, 'cone type');

coneNames = {'L', 'M', 'S'};
for ii = 2 : 4
    fprintf('photoisomerization %s: %d R/sec\n', coneNames{ii-1}, ...
            median(photons(coneType==ii)));
end