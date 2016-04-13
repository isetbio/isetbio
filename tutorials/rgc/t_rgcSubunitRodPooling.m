% t_rgcSubunitRodPooling
% 
% Demonstrates the inner retina object calculation for the subunit RGC
% model (from Gollisch & Meister, 2008, Science; Golisch & Meister, 2010,
% Neuron). 
% 
% Figure 2A of the 2010 Neuron paper shows a simple model of spatial
% pooling for detection sensitivity.  We implement that model here, which
% is only isomerizations, temporal filtering, half wave rectification,
% summation (bipolar).
%
% We explore the properties of that model for different receptor models,
% parameters of the models, simple stimuli, and so forth.  The point of
% this is to get us into the mode of implementing models in the literature
% to try to replicate what is in the papers.
% 
% 3/2016 BW JRG HJ (c) isetbio team

%%
ieInit

%% Movie of the a monochromatic region of rod absorptions (TODO)


% For the moment, we take what's up there which is a bunch of cone
% isomerizations.  We will make a directory for creating rod stimuli as
% well.

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

% Grating subunit stimulus
% params.barWidth = 24;
% iStim = ieStimulusGratingSubunit;
% absorptions = iStim.absorptions;

% White noise
% iStim = ieStimulusWhiteNoise;

% Show raw stimulus for osIdentity
coneImageActivity(absorptions,'dFlag',true);

%% Photocurrent calculation

% We have the moment by moment absorptions.  We now want to create the
% tempmorally filtered version.  If G&M had given us a temporal impulse
% response for the photoreceptor, we would have used it.  For this
% calculation, we use the ISETBIO default.
os = osCreate('linear');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
os = osSet(os, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
os = osSet(os, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
os = osCompute(os,absorptions);

% Plot the photocurrent for a pixel
% osPlot(os,absorptions);
% 
% Can we make a movie of the photocurrent over time?
%

%% G&M tell us to half-wave rectify the photocurrent response.

% What the hell does that mean?  Where is zero?
% We should write a hwRect function
%
%   out = hwRect(data,val);
%
val = -50;
hwrCurrent = max(os.coneCurrentSignal,val);


%% Then we need a little spatial summation over the cones

% This is like a bipolar cell, but actually it could be the same code as in
% the spatial summation of the RGC
%
%  out = spatialSummation(os.coneCurrentsignal,params);
%
%  spatialTemporalSummation()
%

osSize = size(hwrCurrent)

% Set subunit size
% When numberSubunits is set to the RF size, every pixel is a subunit
% This is the default, after Gollisch & Meister, 2008
% sRFcenter = mosaicGet(innerRetina0.mosaic{1},'sRFcenter');
% mosaicSet(innerRetina0.mosaic{1},'numberSubunits',size(sRFcenter));

% Alternatively, have 2x2 subunits for each RGC
% mosaicSet(innerRetina0.mosaic{1},'numberSubunits',[2 2]);

numberSubunits = [2 2];%mosaic.numberSubunits;
suSize1 = floor(osSize(1)/numberSubunits(1));
suSize2 = floor(osSize(2)/numberSubunits(2));
suCtr = 0;
for suInd1 = 1:numberSubunits(1)
    for suInd2 = 1:numberSubunits(2)
        suCtr = suCtr+1;
        xsm = (suInd1-1)*suSize1 + 1: (suInd1)*suSize1;
        ysm = (suInd2-1)*suSize2 + 1: (suInd2)*suSize2;
        
        subunitResponseTemp = hwrCurrent(xsm,ysm,:);
        subunitResponseRS = reshape(subunitResponseTemp,[length(xsm)*length(ysm),osSize(3)]);
        fullResponseSmall(suCtr,:) = mean(subunitResponseRS,1);
    end
end

% fullResponse{xcell,ycell,1} = mean(mosaic.rectifyFunction(fullResponseSmall));






%% Build the inner retina object

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(os, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','subunit','type','on midget');
innerRetina0.mosaicCreate('model','lnp','type','on midget');

irPlot(innerRetina0,'mosaic');

% Set subunit size
% When numberSubunits is set to the RF size, every pixel is a subunit
% This is the default, after Gollisch & Meister, 2008
sRFcenter = mosaicGet(innerRetina0.mosaic{1},'sRFcenter');
mosaicSet(innerRetina0.mosaic{1},'numberSubunits',size(sRFcenter));

% Alternatively, have 2x2 subunits for each RGC
% mosaicSet(innerRetina0.mosaic{1},'numberSubunits',[2 2]);
%% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, os);
irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Show me the PSTH for one particular cell

% irPlot(innerRetina0, 'psth response','cell',[2 2]);
% irPlot(innerRetina0, 'raster','cell',[1 1]);