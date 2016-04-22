% t_rgcSubunit
% 
% Demonstrates the inner retina object calculation for the subunit RGC
% model (from Gollisch & Meister, 2008, Science).
% 
% This is a simplistic implementation of a bipolar-like subunit model for
% RGC computation. The receptive field is broken up into a number of
% subunit fields; at each time step, the input to each subunit is
% summed linearly, and the subunits activations are half-wave rectified and
% summed. The original Gollisch & Meister model is meant to account for
% latencies of spikes after a grating presentation, and the implementation
% here attaches the subunit model as a front end to the spike generating
% code by Pillow et al., Nature, 2008.
% 
% 3/2016 BW JRG HJ (c) isetbio team

%%
clear
ieInit

%% Movie of the cone absorptions 
% % Get data from isetbio archiva server
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% a = rd.listArtifacts;
% 
% % Pull out .mat data from artifact
% whichA =1 ;
% data = rd.readArtifact(a(whichA).artifactId);
% % iStim stores the scene, oi and cone absorptions
% iStim = data.iStim;
% absorptions = iStim.absorptions;

%% Grating subunit stimulus
clear params
params.fov = 1; % degrees
params.barWidth = 1;
iStim = ieStimulusGratingSubunit(params);
absorptions = iStim.absorptions;

%% White noise
% iStim = ieStimulusWhiteNoise;

%% Show raw stimulus for osIdentity
figure;
for frame1 = 1:size(iStim.sceneRGB,3)
    imagesc(squeeze(iStim.sceneRGB(:,:,frame1,:)));
    colormap gray; drawnow;
end
close;

% coneImageActivity(iStim.absorptions,'dFlag',true);

%% Outer segment calculation
% 
% Input = RGB
osLinear = osCreate('linear');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osLinear = osSet(osLinear, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osLinear = osSet(osLinear, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
% osI = osSet(osI, 'rgbData', 0.6-iStim.sceneRGB);
osLinear = osCompute(osLinear,absorptions);

% % Plot the photocurrent for a pixel
osPlot(osLinear,absorptions);

%% Plot cones
cone_mosaic = absorptions.human.coneType;
[xg yg] = meshgrid([1:123,1:150]);
xg2 = xg(1:123,1:150); yg2 = yg(1:123,1:150);

figure; scatter(xg2(:),yg2(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
%         
%% Build the inner retina object

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina1 = irCreate(osLinear, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina1.mosaicCreate('model','lnp','type','off parasol');
irPlot(innerRetina1,'mosaic');

params.eyeRadius = 3;        % Radius in mm
innerRetina2 = irCreate(osLinear, params);
innerRetina2.mosaicCreate('model','lnp','type','off midget');
hold on;
irPlot(innerRetina2,'mosaic');

%% Build the bipolar cell filter

% osFilter = osGet(osI,'mConeFilter');
osFilter = osLinear.mConeFilter;

osFilterDerivativeShort = diff(osFilter);
osFilterDerivative = interp1(1:length(osFilterDerivativeShort),osFilterDerivativeShort,1:length(osFilter))';

bipolarFilter = osFilter + .5*osFilterDerivative;

figure; 
plot(bipolarFilter);
hold on;
plot(osFilter,'g');
plot(osFilterDerivative,'r');
xlabel('Time (sec)'); ylabel('pA / (R*/sec)');
set(gca,'fontsize',16);

legend('bipolar f(t)','os f(t)','os f''(t)');

%% Filter to get bipolar cell output

eZero = -50;
hwrCurrent = ieHwrect(osLinear.coneCurrentSignal,eZero);

kernel = fspecial('gaussian',[9,9],3);

bipolar = ieSpaceTimeFilter(hwrCurrent,kernel);

strideSubsample = 4;
bipolarSubsample = ieImageSubsample(bipolar, strideSubsample);

szBS = size(bipolarSubsample);
bipolarSubsampleRS = reshape(bipolarSubsample,szBS(1)*szBS(2),szBS(3));

bipolarOutputRS = convn(bipolarSubsampleRS,bipolarFilter','full');

bipolarOutput = reshape(bipolarOutputRS,szBS(1),szBS(2),size(bipolarOutputRS,2));


%% Create outer segment identity in order to pass bipolarOutput to innerRetina

% Input = RGB
osI = osCreate('displayRGB');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
bipolarOutputRGB = repmat(bipolarOutput./3,[1 1 1 3]);
osI = osSet(osI, 'rgbData', bipolarOutputRGB);
%% Build the inner retina object

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(osI, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','lnp','type','off parasol');

irPlot(innerRetina0,'mosaic');

% Set subunit size
% When numberSubunits is set to the RF size, every pixel is a subunit
% This is the default, after Gollisch & Meister, 2008
sRFcenter = mosaicGet(innerRetina0.mosaic{1},'sRFcenter');
% mosaicSet(innerRetina0.mosaic{1},'numberSubunits',size(sRFcenter));

% Alternatively, have 2x2 subunits for each RGC
% mosaicSet(innerRetina0.mosaic{1},'numberSubunits',[2 2]);
%% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, osI);
irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Show me the PSTH for one particular cell

% irPlot(innerRetina0, 'psth response','cell',[2 2]);
% irPlot(innerRetina0, 'raster','cell',[1 1]);