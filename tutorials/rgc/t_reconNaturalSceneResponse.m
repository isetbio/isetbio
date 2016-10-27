% t_reconNaturalSceneResponse

% This tutorial generates an RGC mosaic response to a natural scene
% stimulus. The response can be used with the repository at 
% github.com/Chichilnisky-Lab/RGC-Reconstruction to generate a
% reconstructed image.
%
% The response is generated using the isetbio framework. The general steps 
% are as follows:
%   * display
%   * scene
%   * optical image [oi]
%   * sensor        [cone absorptions]
%   * outer segment [cone current]
%   * inner retina  [rgc mosaics]
% 
% The scene is set to be a natural image displayed on a particular
% calibrated monitor.  The standard RGC LNP and GLM models predict responses
% based on the stimulus RGB values, so both the optical image structure
% that computes the image that falls on the retina and the sensor structure 
% that computes cone absorptions are not used in the computation of the RGC 
% response. The  oi and sensor are created to generate a cone mosaic
% of the appropriate size in order to build an RGC. A special type of outer
% segment is created to store the RGB values instead of the cone photon
% absorptions. The RGC responses are then computed from the outer segment
% object RGB values.

%% Parameters

natScenStim = 1; % if 0, execute moving bar

blankOnsetFrames = 20;
stimulusFrames = 10;

%% Display
display = displayCreate('LCD-Apple');

%% Scene
if ~natScenStim
    % moving bar stimulus, show spatial structure in response
    scene = sceneCreate;
    % scene = sceneCreate('rings rays');
    params.barWidth = 10;
else
    % INSERT NATURAL IMAGE HERE
    % Need to resize
    Ibig = imread('peppers.png');
    I = rgb2gray(imresize(Ibig,0.5));
    
    params.meanLuminance = 200;
    scene = sceneFromFile(I, 'rgb', params.meanLuminance, display);
end
%% Optical image

oi  = oiCreate('wvf human');

frameTotal = blankOnsetFrames + stimulusFrames;

params.fov = 1.5;
% params.barWidth = 10;

sceneRGB = zeros([sceneGet(scene, 'size'), frameTotal, 3]);

%% Cone mosaic

sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, params.fov, scene, oi);

params.expTime = (1/125);
params.timeInterval = (1/125);%  
sensor = sensorSet(sensor, 'exp time', params.expTime); 
sensor = sensorSet(sensor, 'time interval', params.timeInterval); 

for frameNumber = 1:frameTotal
    
    if ~natScenStim
    barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;
    
    % Bar at this time
        if frameNumber > blankOnsetFrames
            colStart = frameNumber + 1 - blankOnsetFrames;
            colEnd   = colStart + params.barWidth - 1;
            % barMovie(:,t-startFrames + 1:(t-startFrames+1+params.barWidth-1),:) = 1;
            barMovie(:,colStart:colEnd,:) = 1;
        end
        scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);
        
    end
    
    %     scene = sceneCreate('rings rays');
    %    scene = sceneFromFile(I, 'rgb');
    
    scene = sceneSet(scene, 'h fov', params.fov);

    % Get scene RGB data
    if frameNumber < blankOnsetFrames || frameNumber > blankOnsetFrames + stimulusFrames
        sceneRGB(:,:,frameNumber,:) = 0.5*ones(size(sceneGet(scene,'rgb')));
    else        
        sceneRGB(:,:,frameNumber,:) = sceneGet(scene,'rgb');
    end
    
    % Compute optical image
    % oi = oiCompute(oi, scene);
    
end

%% Outer segment 

os = osCreate('displayrgb');

sceneSize = sceneGet(scene,'size');
retinalPatchWidth = sensorGet(sensor,'width','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

os = osSet(os, 'patchSize', retinalPatchWidth);

timeStep = (1/125)/1;
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', sceneRGB);

retinalPatchSize = osGet(os,'size');

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 4;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetina = irCreate(os,paramsIR);
innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
innerRetina = irCompute(innerRetina,os);

%% Plot outputs
% figure; plot(RGB2XWFormat(innerRetina.mosaic{1}.responseLinear));
psth = innerRetina.mosaic{2}.get('psth');
figure; plot(RGB2XWFormat(psth)');
figure; ieMovie(psth(:,:,1:100:end));