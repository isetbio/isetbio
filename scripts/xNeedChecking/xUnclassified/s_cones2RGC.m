%% s_rgcCones2RGC
%
% The calculation starting from cone array data to RGC spikes.
%
%
% See also:  s_rgcScene2Cones
%
% (c) Stanford VISTA Team

%% First time through, at least, initialize ISET.  Not always necessary
s_initISET; 

%% Load a cone image pre-computed by s_rgcScene2Cones
%fName = 'edgeEyeMvmnts';
fName = 'eagleEyeMvmnts';
load(fullfile(isetbioRootPath,'tmp',fName));

% To visualize do this
%   scene = cones.scene;
%   vcAddAndSelectObject(scene); sceneWindow;

%% Create RGC structure

% creating the RGC parameter object
rgcP = rgcParameters;

% The absorptions structure is also the RGC parameter data
rgcP.set('scene',cones.scene);
rgcP.set('oi',cones.oi);
rgcP.set('sensor',cones.sensor);
rgcP.set('cone voltages',cones.data);

% What is the default?  RGC spacing the same as cone spacing (I think).
rgcP.addLayer('on parasol', 20);  
rgcP.addLayer('off parasol');  
 
rgcP.addLayer('on midget');  
rgcP.addLayer('off midget');  

rgcP.addLayer('small bistratified');  % Sign of signal seems wrong. Check.

% Visualize some stuff
%  whichLayer = 1;
%  rgcVisualize('all tirf',rgcP, whichLayer)
%  rgcVisualize('fbtr',rgcP, whichLayer)
%  rgcVisualize('cptr',rgcP, whichLayer)
%  rgcVisualize('RF Mesh',rgcP)
%
%  rgcVisualize('RF image',rgcP, whichLayer);  % Single RGC receptive field 
%  rgcVisualize('RF mesh',rgcP, whichLayer);   % Single RGC receptive field 
%  rgcVisualize('RF mesh',rgcP, 1);   % Single RGC receptive field 
%  rgcVisualize('RF mesh',rgcP, 2);   % Single RGC receptive field 
%
% Would be nice to show the 1 SD outline over the cone mosaic.
% How about it guys?
%
%  rgcVisualize('RF center',rgcP, whichLayer);         % Just the center
%  rgcVisualize('RF surround',rgcP, whichLayer);       % Just the surround
%
%  rgcVisualize('tirf',rgcP, whichLayer);     % Temporal Impulse Response

%% RGC Spike computation

 
% If we don't reinitialize rgcP, we get nothing.  So check the recompute
% aspects of this.  Probably we should also clear out earlier data ...

% TODO TO THINK
% OOOPS.  There is a problem with the noise calculation.
%
% When we return, we want to understand how the voltage levels for noise
% and threshold are influencing the behavior of the spiking outputs. We
% also want to be able to try different models for how the spike voltages
% add together (see rgcComputeSpikes).
%
% We want to move on to the SVM classification stuff.
%
% Change the parameter names in sets/gets to be better.  Add coneVolts, get
% rid of the 'current' part, eliminate weird names in the mapping that we
% would never use.  
% We want to get the spike rate proper by figuring out the entire stimulus
% duration and add that to the get part of the layer.

nL = rgcP.get('nL');
layer = cell(1,nL);
for ii=1:nL, layer{ii} = rgcP.get('layer',ii); end

hasFeedback = ones(1,nL);
for ii=1:nL, layer{ii}.hasFeedback = hasFeedback(ii); end

hasCoupling = zeros(1,nL);
for ii=1:nL, layer{ii}.hasCoupling = hasCoupling(ii); end

% Spike threshold as a percentage of the RGC voltage swing
vSwingOriginal = layer{1}.get('vswing');
for ii=1:nL, layer{ii}.set('rgc volt thresh',0.2*vSwingOriginal); end

% layer{1}.set('cone weights',[.3 0 0; 0 0 .3]);
% layer{1}.get('cone weights')

% We also want the largest value in the coupling and feedback to be smaller
% than the spike threshold.  So attend to that around here.
% fbTimeOriginal = layer{1}.get('fbtr');
% for ii=1:nL
%     layer{ii}.set('rgc volt thresh',0.5*fbTimeOriginal);
% end

rgcComputeSpikes(rgcP);

% Things have changed in rgcP, too.  You can see them using:
%   rgcP.layers{1}.get('rgc volt thresh')
%   rgcP.layers{1}.get('range lin ts')
%

%% **** Movies (time series) ********

% question: are absorptions stored somewhere and are they different from volts?
% rgcVisualize('Cone Absorptions',rgcP);  % NYI. Cones absorption movie. 

% rgcVisualize('Cone Voltages',rgcP);     % Cones voltage movie

%
% [rgcComputeSpikes]: Layer 1 (on parasol)
% [rgcComputeSpikes]: Layer 2 (off parasol)
% [rgcComputeSpikes]: Layer 3 (on midget)
% [rgcComputeSpikes]: Layer 4 (off midget)
% [rgcComputeSpikes]: Layer 5 (small bistratified)

% Spikes
for whichLayer = 3:4
    % rgcVisualize('component TS',rgcP, whichLayer);    % RGC, movie of separate center and surround linear time series
    % rgcVisualize('Lin TS',rgcP, whichLayer)             % RGC linear t-series movie (driving current)
    % rgcVisualize('RGC TS',rgcP, whichLayer)             % RGC t-series movie (driving current + feedback + coupling)
    rgcVisualize('spikes',rgcP, whichLayer);            % RGC, spiking movie
end

%% RGC time series
for whichLayer = 3:4
    % rgcVisualize('component TS',rgcP, whichLayer);    % RGC, movie of separate center and surround linear time series
    % rgcVisualize('Lin TS',rgcP, whichLayer)             % RGC linear t-series movie (driving current)
    rgcVisualize('RGC TS',rgcP, whichLayer)             % RGC t-series movie (driving current + feedback + coupling)
    % rgcVisualize('spikes',rgcP, whichLayer);            % RGC, spiking movie
end


%% **** Images (TS mean) ***********
% rgcVisualize('Cone Volts mean',rgcP);   % Cones voltage image  (mean of ts)
% 

for whichLayer = 1:nL
    rgcVisualize('Lin TS mean',rgcP, whichLayer)      % RGC linear t-series(driving current), (mean of ts)    
    %rgcVisualize('RGC TS mean',rgcP, whichLayer)     % RGC t-series membrane voltage (driving current + feedback + coupling), (mean of ts)    
    %rgcVisualize('Spikes TS mean',rgcP, whichLayer); % RGC, spike rate, (mean of ts)
end

%% Feedback and coupling

for whichLayer = 1:nL
    rgcVisualize('all tirf',rgcP, whichLayer)        % RGC linear temporal response function   
    rgcVisualize('feedback tr',rgcP, whichLayer)     % RGC feedback temporal response function   
    rgcVisualize('coupling tr',rgcP, whichLayer)     % RGC coupling temporal response function   
end


%% Get some parameters
spkTS = cell(1, nL);
for whichLayer = 1:nL
    spkTS{whichLayer} = layer{whichLayer}.get('spike time series');
end
%% Notes: absorptions is a structure with a lot of fields
%
% The data field is volts that are adapated
%     unadapated: volts prior to adaptation
%     adaptGain:  Scale factor applied to unadapated to get data
%     sensor:     The volts field appears to have the noise-free mean
%     oi and scene - Bundled for later RGC calculations?
%
% The adaptGain contain the scale factors for the different cone classes.
% We always make the gain for the empty spaces ('k' filter) 1.  
% These K values should probably be removed or all set to zero and thus
% meaningless at some point.absorptions.adaptGain(:)

% Visualize the mean adapted and unadapted data
%
% sensorAdapted = sensor; 
%
% Adaptation is turned off ...
%  sensorUnadapted = sensorSet(sensorAdapted,'name','No adaptation');
%  sensorUnadapted = sensorSet(sensorUnadapted,'volts',mean(absorptions.unadapted,3));
%  vcAddAndSelectObject(sensorUnadapted); sensorImageWindow;
%
% Adaptation is turned on
%  sensorAdapted = sensorSet(sensorAdapted,'name','adaptation');
%  sensorAdapted = sensorSet(sensorAdapted,'volts',mean(absorptions.data,3));
%  vcAddAndSelectObject(sensorAdapted); sensorImageWindow;

% This is the mean number of absorptions in a frame
%  figure; imagesc(mean(absorptions.data,3)); colormap(gray); axis image
%  fprintf('Time samples: %0.4f s\n',sensorGet(sensor,'expTime'))
%  colorbar

% End
