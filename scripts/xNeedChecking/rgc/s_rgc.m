
clear

% % % load('sensor_short_CV_3_Cont_6000_fov_08.mat')

%%  preliminaries - define gabor color opponent stimulus, scene, oi, sensor

% build scene
params.image_size = 64;
params.meanLuminance = 100;
params.nsteps = 10;
params.fov = 0.8;
[scene, display] = sceneHorwitzHassWhiteNoise(params);
% displayClose;
% % build optical image
oi  = oiCreate('wvf human');

% build sensor for white noise
sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);

%%
% build outersegment
identityOS = osCreate('identity');

fprintf('Computing RGB scene data    ');
for step = 1:params.nsteps
    fprintf('     \n');
    fprintf('\b\b\b%02d%%', round(100*step/params.nsteps));
    fprintf('     \n');
    [scene, display] = sceneHorwitzHassWhiteNoise(params);
    sceneRGB(:,:,:,step) = sceneGet(scene,'rgb');
    
end
%%
identityOS = osSet(identityOS, 'rgbData', sceneRGB);

%%
% build rgc

% rgc1 = rgcLinear(sensor, osIdentity, 'right', 3.75, 180);
% rgc1 = rgcLNP(sensor, osIdentity, 'right', 3.75, 180);

rgc1 = rgcGLM(sensor, osIdentity, 'right', 3.75, 180);

%%

% % % tic
rgc1 = rgcCompute(rgc1, sensor, identityOS);
% % % toc
%
% % % tic
rgcPlot(rgc1, sensor, identityOS);
% % % toc
% 
% %%
% 
rgcMovie(rgc1, identityOS);