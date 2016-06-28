% t_rgcOsComparison
% 
% Demonstrates the inner retina object in isetbio. For a stimulus that
% consists of a bar that sweeps from left to right, the response of an rgc
% mosaic is computed based on the raw RGB data using the osIdentity outer
% segment object, the linear outer segment current response using osLinear,
% and the biophysical outer segment current response using osBioPhys.
% 
% 3/2016 BW JRG HJ (c) isetbio team

%%
ieInit

%% Movie of the cone absorptions 
% Get data from isetbio archiva server
rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
a = rd.listArtifacts;

% Pull out .mat data from artifact
whichA =1 ;
data = rd.readArtifact(a(whichA).artifactId);
% iStim stores the scene, oi and cone absorptions
iStim = data.iStim;
absorptions = iStim.absorptions;


%% Show raw stimulus for osIdentity
figure;
for frame1 = 1:size(iStim.sceneRGB,3)
    imagesc(squeeze(iStim.sceneRGB(:,:,frame1,:)));
    colormap gray; drawnow;
end
close;

%% Outer segment calculation
% 
% Input = RGB
osI = osCreate('identity');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
osI = osSet(osI, 'rgbData', iStim.sceneRGB);
% os = osCompute(sensor);

% % Plot the photocurrent for a pixel
% osPlot(osI,absorptions);
%% Build the inner retina object

params.name = 'macaque phys'
params.outersegment = osI;
params.eyeSide = 'left'; 
params.eyeRadius = 9; 
params.eyeAngle = 90;
% rgc2 = rgcCreate('rgcPhys', params);
innerRetina0 = irPhys(osI, params);
innerRetina0 = irSet(innerRetina0,'numberTrials',20);

% irPlot(innerRetina0,'mosaic');

%% Compute RGC mosaic responses

% innerRetina0 = irCompute(innerRetina0, osI);

innerRetina0 = irComputeContinuous(innerRetina0, osI);

nTrials = 20;
for tr = 1:nTrials
    innerRetina0 = irComputeSpikes(innerRetina0, osI);
end

psth = mosaicGet(innerRetina0.mosaic{1},'responsePsth');

irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'voltage');
% irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Show me the PSTH for one particular cell

irPlot(innerRetina2, 'psth response','cell',[2 2]);
title('OS Biophys and Coupled GLM');

%%

irPlot(innerRetina2, 'raster','cell',[1 1]);