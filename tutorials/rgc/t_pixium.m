% t_pixium 
% 
% A first-pass attempt at simulating the activations of retinal
% ganglion cells with an array of stimulating electrodes. The stimulus
% image is defined, the electrode array is generated, the electrode
% activations are computed, the RGC mosaics are generated, the RGC mosaic
% responses are computed and the stimulus is inferred from the RGC mosaic
% resposnes using simple linear summation of the STAs.
% 
% 3/2016 JRG (c) isetbio

% % % % % % % % % 
% folks,
% 
% My summary and suggestion about next steps regarding phosphene modeling.
% 
% Our goal is to understand whether a phosphene model based on primate
% retina can be valuable. We want to have a sense of this before investing
% huge effort.
% 
% From Georges, we know that temporal properties of retinal activation are
% very tricky and are the subject of study in the Palanker lab now. This
% will take time. I suggest the following in the meantime. A simplistic
% phosphene simulator consists of two parts:
% 
% 1) The forward transformation array->RGCs is a spatial spread function
% (Gaussian) for each electrode, that randomly covers an interspersed
% collection of cell bodies of different RGC types, and activates each one
% with a sigmoidal probability as a function of current from the electrode.
% The spatial spread of every electrode, and the magnitude of its output,
% vary from one electrode to the next. From Georges we will get a rough
% number for the inter-electrode variability of this spatial spread and
% magnitude. This is sufficient to produce an estimated mapping
% array->RGCs.
% 
% 2) The reconstruction RGCs->perception is done as follows. James uses the
% RF models he has in ISETBIO to do a simplistic static image
% reconstruction (no dynamics) from RGC spikes. This is done by placing a
% copy of the STA in the reconstruction for every RGC spike.  Again, no
% temporal stuff at all.
% 
% From the above, we have array->RGCs and RGCs->perception.   Finally, let
% us assume that we take an image and activate the array in proportion to
% image intensity (in other words, no pre-processing of the image, so
% image->array is the identity transformation).  With the above items, we
% have image->perception.
% 
% I think James and Cordelia should work together on this and start looking
% at the resulting images, varying these parameters:
% 
% - inter-electrode variability of spatial spread/magnitude of current -
% RGC density (i.e. eccentricity) - relative sensitivity of different RGC
% types to electrical stimulation
% 
% This should be pretty easy and seems to me the first tool we would need
% to build intuitions about potential for spatial processing. It will not
% be hard to build temporal properties on top of this. We will save the big
% guns on reconstruction (Nishal) for when we outgrow this first-pass
% approach.
% 
% ej

% 1.5 orders of magnitude variance in threshold of RGC activation

%% Initialize
clear;
ieInit;

% Set the size of implant pixels
electrodeArray.width = 70e-6; % meters
% electrodeArray.width = 140-6; % meters

%% Load image

% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = 1;
params.row = 100;
params.col = 100;
params.fov = 0.7;
% params.vfov = 0.7;
movingBar = ieStimulusBar(params);

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('identity');

retinalPatchWidth = sensorGet(movingBar.sensor,'width','m');
% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

timeStep = sensorGet(movingBar.sensor,'time interval','sec');
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', movingBar.sceneRGB);
% os = osCompute(sensor);

% % Plot the photocurrent for a pixel
% osPlot(os,sensor);

retinalPatchSize = osGet(os,'size');
numberElectrodes = retinalPatchWidth/electrodeArray.width;
%% Build electrode array
% Define the electrode array structure/object

% Size stores the size of the array of electrodes
electrodeArray.size = [numberElectrodes numberElectrodes];

% Builds the matrix of center coordinates for each electrode
% electrodeArray.center(xPos,yPos,:) = [xCoord yCoord];
for xPos = 1:numberElectrodes
    for yPos = 1:numberElectrodes
        electrodeArray.center(xPos,yPos,:) = [(retinalPatchWidth/2)*(xPos-1) + retinalPatchWidth, (retinalPatchWidth/2)*(yPos-1) + retinalPatchWidth];
    end
end

% Build the current stimulation activation window
% Gaussian activation from center of electrode
activationWindow = round(retinalPatchSize(2)/numberElectrodes);
electrodeArray.activation = fspecial('Gaussian', activationWindow, activationWindow/8);

% Visualize Gaussian activation
% figure; imagesc(electrodeArray.activation); 
figure; surf(electrodeArray.activation); 
xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);
%% Compute electrode activations from image

% Get the full image/movie from the identity outersegment
fullStimulus = osGet(os,'rgbData');

% Find electrode activations by taking mean within window
for frame = 1:params.nSteps
    for xPos = 1:numberElectrodes
        for yPos = 1:numberElectrodes
            % Xcoords of window for stimulus
            imageCoordX1 = (activationWindow)*(xPos-1)+1;
            imageCoordX2 = (activationWindow)*(xPos);
            
            % Ycoords of window for stimulus
            imageCoordY1 = (activationWindow)*(yPos-1)+1;
            imageCoordY2 = (activationWindow)*(yPos);
            
            % Pull out piece of stimulus and take mean
            electrodeStimulus = squeeze(fullStimulus(imageCoordX1:imageCoordX2,imageCoordY1:imageCoordY2,frame,:));
            electrodeArray.activation(xPos,yPos,frame) = mean(electrodeStimulus(:));
        end
    end
end

% Apply Gaussian

%% Build RGC array

clear params
params.name    = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 12;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

model   = 'pool';    % Computational model
innerRetina = irCreate(os,params);
innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);

innerRetina = irCompute(innerRetina, os);

%% Build RGC activation functions

%% Compute RGC activations

%% Invert representation to form image/movie

