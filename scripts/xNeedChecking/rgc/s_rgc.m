
% s_rgc
% 
% This script builds scene, sensor, outersegment and rgc objects in isetbio
% and computes the RGC responses. The RGC responses are stored in the rgc
% mosaic object and can be viewed as a movie by running the command on the
% last line.
% 
% isetbio
% 9/2015 JRG

clear

%% build scene
params.image_size = 96;
params.meanLuminance = 100;
params.nsteps = 30;
params.fov = 0.6;
[scene, display] = sceneHorwitzHassWhiteNoise(params);

% vcAddObject(scene); sceneWindow

% displayClose;

% % build optical image
oi  = oiCreate('wvf human');

%% build sensor for white noise
params.nsteps = 1;
sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
params.nsteps = 30;
%% build outersegment
identityOS = osCreate('identity');

% sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);
sceneRGB = sceneHorwitzHassBarRGB(params);
% sceneRGB = 0.5*ones(params.image_size, params.image_size, params.nsteps, 3);

% for frame = 1:params.nsteps
%     params.freq =  [5 ]; % spatial frequencies of 1 and 5
%     params.contrast = [0.6]; % contrast of the two frequencies
%     params.ang  = [0,]; % orientations
%     params.ph  = [0 ]; % phase
%     sharmonic = sceneCreate('harmonic',params);
%     sceneRGB(:,:,frame,:) = sceneGet(sharmonic,'rgb');
% end
 
identityOS = osSet(identityOS, 'rgbData', sceneRGB);

%% build rgc

% rgc1 = rgcLinear(scene, sensor, osIdentity, 'right', 3.75, 180);
% rgc1 = rgcLNP(scene, sensor, osIdentity, 'right', 3.75, 180);
% rgc1 = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
% rgc1 = rgcSubunit(scene, sensor, identityOS, 'right', 3.0, 180);

rgc1 = rgcCreate('LNP', scene, sensor, identityOS, 'right', 3.0, 180);
% rgc1 = rgcCreate('Linear', scene, sensor, identityOS, 'right', 3.0, 180);
% rgc1 = rgcCreate('GLM', scene, sensor, identityOS, 'right', 3.0, 180);
% rgcPlot(rgc1, 'mosaic');
%% compute rgc

% % % tic
rgc1 = rgcCompute(rgc1, identityOS);
% % toc

% % tic
rgcPlot(rgc1, 'psthResponse');
% rgcPlot(rgc1, 'rasterResponse');
% % % toc
%% With linear cone response

% linearOS = osCreate('linear');
% linearOS = osCompute(linearOS, sensor);
% % 
% % rgc1 = rgcGLM(scene, sensor, linearOS, 'right', 3.0, 180);
% rgc2 = rgcCreate('glm',scene, sensor, linearOS, 'right', 3.0, 180);
% % 
% rgc2 = rgcCompute(rgc2, linearOS);
% % 
% 
% rgcPlot(rgc2, 'linearResponse');
% rgcPlot(rgc2, 'spikeResponse');
%% build movie
% 
% rgcMovie(rgc1, identityOS);