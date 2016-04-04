% t_pixium 
% 
% A first-pass attempt at simulating the activations of retinal
% ganglion cells with an array of stimulating electrodes. The stimulus
% image is defined, the electrode array is generated, the electrode
% activations are computed, the RGC mosaics are generated, the RGC mosaic
% responses are computed and the stimulus is inferred from the RGC mosaic
% resposnes using simple linear summation of the STAs.
% 
% Outline of computation:
% 1. Load image/movie
% - Moving bar stimulus sweeping left to right
% 
% 2. Outer segment representation
% - The osIdentity object does no computation and stores the stimulus RGB
% data.
%  Check iestimulusbar for size change
% 
% 3. Build electrode array
% - The electrode array is simulated as a square array of a certain size
% with Gaussian weights on their activations. The activation of each
% electrode is the mean of the CHANGE TO HEX
% 
% 4. Compute electrode activations from image/movie
% 5. Build RGC array
% 6. Calculate RGC input - only one, change to input from multiple
% 7. Build RGC activation functions
% 8. Compute RGC activations/spikes
% 9. Invert representation to form image/movie
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
% 
% 1.5 orders of magnitude variance in threshold of RGC activation
% 
% Add 5 Hz stimulation pulses


%% Parameters to alter

% Electrode size
% Retinal eccentricity
% Field of view/stimulus size
% Stimulus length
% Activation curve
% Spatial activation of electrode
% Electrode PWM

%% Initialize
clear;
ieInit;
 
% Set the size of implant pixels
electrodeArray.width = 70e-6; % meters
% electrodeArray.width = 140-6; % meters

% Set horizontal field of view
params.fov = 2.7;

%% Load image

% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = 90;
params.row = 100;
params.col = 100;
% params.fov = 2.7;
% params.vfov = 0.7;
movingBar = ieStimulusBar(params);

%% Show raw stimulus for osIdentity
figure;
for frame1 = 1:size(movingBar.sceneRGB,3)
    imagesc(squeeze(movingBar.sceneRGB(:,:,frame1,:)));
    colormap gray; drawnow;
end
close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('identity');

sceneSize = sceneGet(movingBar.scene,'size');
retinalPatchWidth = sensorGet(movingBar.absorptions,'width','m');
% retinalPatchHeight = sensorGet(movingBar.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

timeStep = sensorGet(movingBar.absorptions,'time interval','sec');
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', movingBar.sceneRGB);
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');
numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+0;
numberElectrodes = numberElectrodesX*numberElectrodesY;
%% Build electrode array
% Define the electrode array structure/object

% Size stores the size of the array of electrodes
electrodeArray.size = [numberElectrodesX numberElectrodesY];

% Build the matrix of center coordinates for each electrode
% electrodeArray.center(xPos,yPos,:) = [xCoord yCoord];
x0 = .5e-4; y0 = .5e-4;
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
        electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/2)*(xPos-1) + electrodeArray.width, y0+(electrodeArray.width/2)*(yPos-1) + electrodeArray.width];
    end
end


th = (0:1/6:1)'*2*pi;
xh = electrodeArray.width/4*cos(th);
yh = electrodeArray.width/4*sin(th);
figure
fill(x,y,'r')
axis square

% % Plot electrode array
eaSize = size(electrodeArray.center);
% figure;
% hold on;
% for i = 1:eaSize(1)
%     for j = 1:eaSize(2)
% %         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
%         plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
%     end
% end

% Build the current stimulation activation window
% Gaussian activation from center of electrode
activationWindow = floor(retinalPatchSize(2)/numberElectrodesX);
electrodeArray.spatialWeight = fspecial('Gaussian', activationWindow, activationWindow/8);

% Visualize Gaussian activation
% figure; imagesc(electrodeArray.spatialWeight); 
figure; surf(electrodeArray.spatialWeight); 
xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);
%% Compute electrode activations from image

% Get the full image/movie from the identity outersegment
fullStimulus = osGet(os,'rgbData');

% Find electrode activations by taking mean within window
for frame = 1:params.nSteps
    for xPos = 1:numberElectrodesX
        for yPos = 1:numberElectrodesY
            % Xcoords of window for stimulus
            imageCoordX1 = (activationWindow)*(xPos-1)+1;
            imageCoordX2 = (activationWindow)*(xPos);
            
            % Ycoords of window for stimulus
            imageCoordY1 = (activationWindow)*(yPos-1)+1;
            imageCoordY2 = (activationWindow)*(yPos);
            
            if imageCoordX2 > size(fullStimulus,2); imageCoordY2 = size(fullStimulus,2); end;
            if imageCoordY2 > size(fullStimulus,1); imageCoordY2 = size(fullStimulus,1); end;
            % Pull out piece of stimulus and take mean
            electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,frame,:));
            electrodeArray.activation(xPos,yPos,frame) = mean(electrodeStimulus(:));
            % imagesc(electrodeStimulus); title(sprintf('%2.2f',mean(electrodeStimulus(:))));
        end
    end
end

% Show plot

% Apply Gaussian

%% Add 5 Hz spiking of stimulus

%% Build RGC array

clear paramsIR
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 12;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetina = irCreate(os,paramsIR);
innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);

irPlot(innerRetina,'mosaic');
% % figure;
% hold on;
% for spInd = 1:length(innerRetina.mosaic)
% for i = 1:eaSize(1)
%     for j = 1:eaSize(2)
%         subplot(floor(sqrt(length(innerRetina.mosaic))),floor(sqrt(length(innerRetina.mosaic))),spInd); 
%         hold on;
% %         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
%         plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
%     end
% end
% end

metersPerPixel = retinalPatchWidth/retinalPatchSize(2);

% innerRetina = irCompute(innerRetina, os);

% Plot the electrode array over spatial receptive fields
% eaSize = size(electrodeArray.center);
% hold on;
% for i = 1:eaSize(1)
%     for j = 1:eaSize(2)
%         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
%     end
% end
%% Calculate RGC input
% Weight electrode activation by Gaussian as a function of distance between
% centers


for frame = 1:params.nSteps
    for mosaicInd = 1:length(innerRetina.mosaic)
        [xc yc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
        for xind = 1:xc
            for yind = 1:yc
                % Find closest electrode
                % Electrode centers
                electrodeCenter = reshape(electrodeArray.center,[numberElectrodes,2]);
                % RGC center
                rgcCenter = innerRetina.mosaic{mosaicInd}.cellLocation{xind,yind}*metersPerPixel;
                % Find distance
                centerDistanceCoords = repmat(rgcCenter,[numberElectrodes,1]) - electrodeCenter;
                % Weight electrode activation by Gaussian according to distancej
                centerDistance = bsxfun(@(x,y)sqrt(x.^2+y.^2),centerDistanceCoords(:,1),centerDistanceCoords(:,2));
                centerDistanceRS = reshape(centerDistance,[numberElectrodesX,numberElectrodesY]);
                
                [minDistance, minDistanceInd] = min(centerDistance);
                [xmin,ymin] = ind2sub([numberElectrodesX,numberElectrodesY],minDistanceInd);
                minXY(yind,xind,frame,mosaicInd,:) = [xmin ymin];
                innerRetinaInput(yind,xind,frame,mosaicInd) = electrodeArray.activation(xmin,ymin,frame)*exp(-minDistance/2e-4);
                
%                 for xind2 = 1:numberElectrodesX
%                     for yind2 = 1:numberElectrodesY
%                         innerRetinaInput(yind,xind,frame,mosaicInd) = ...
%                             innerRetinaInput(yind,xind,frame,mosaicInd) + ...
%                             electrodeArray.activation(xind2,yind2,frame)*exp(-centerDistanceRS(xind2,yind2)/2e-4);
%                     end
%                 end
%                 
            end
        end
    end
end


% figure; imagesc(squeeze(minXY(:,:,10,1,1)))
% figure; imagesc(squeeze(minXY(:,:,10,1,2)))
%% Build RGC activation functions

% figure; hold on;
% for i = 1:10
% x0 = -.5*rand(1,1); thr = 10*rand(1,1);
% % figure; 
% xp = -1:.01:3;
% plot(-x0+(xp),1./(1+exp(-thr*(x0+(xp)))));
% end

% figure; hold on;
for mosaicInd = 1:length(innerRetina.mosaic)
    [xc yc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
    for xind = 1:xc
        for yind = 1:yc
            thr = 20*rand(1,1);
            i0 = -.5*rand(1,1);
            % innerRetinaFunction{xind,yind,mosaicInd}= @(iElectrode) 1./(1+exp(-thr*(i0+iElectrode)));
            innerRetinaThreshold{xind,yind,mosaicInd} = [thr i0];
            % plot( innerRetinaFunction{xind,yind,mosaicInd}(-1:.1:4));
        end
    end
end

% Plot something
%% Compute RGC activations
for mosaicInd = 1:length(innerRetina.mosaic)
    clear innerRetinaActivation 
    [xc yc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
    for xind = 1:xc
        for yind = 1:yc
            
            for frame = 1:params.nSteps
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction{xind,yind,mosaicInd}(innerRetinaInput(xind,yind,mosaicInd));
                funcParams = innerRetinaThreshold{xind,yind,mosaicInd};
                thr = funcParams(1); i0 = funcParams(2);
                % innerRetinaFunction = @(iElectrode) 10./(1+exp(-thr*(i0+iElectrode)));
                % innerRetinaFunction = @(iElectrode) log(1./(1+exp(-thr*(i0+iElectrode))));
                innerRetinaFunction = @(iElectrode) (5*iElectrode);
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction(innerRetinaInput(xind,yind,mosaicInd));
                
                innerRetinaActivation{xind,yind}(frame) = 10*innerRetinaFunction(innerRetinaInput(xind,yind,frame,mosaicInd));
            end
            mosaicSet(innerRetina.mosaic{mosaicInd},'responseLinear', innerRetinaActivation);
            
        end
    end
end

%% Compute RGC spiking
numberTrials = 10;
for tr = 1:numberTrials
    innerRetina = irComputeSpikes(innerRetina);
end
%% Invert representation to form image/movie

irReconstruct(innerRetina);
