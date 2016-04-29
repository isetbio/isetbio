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
% 
% 1.5 orders of magnitude variance in threshold of RGC activation
% 
% Add 5 Hz stimulation pulses




%% Initialize
clear;
ieInit;

%% Parameters to alter

% Electrode size
% Set the size of implant pixels
% electrodeArray.width = 70e-6; % meters
electrodeArray.width = 140e-6; % meters

% Retinal patch eccentricity
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 2.7;

% Stimulus length
nSteps = 80;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 


%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 100;
params.col = 100;
params.fov = fov;
% % params.vfov = 0.7;
% movingBar = ieStimulusBar(params);


%%% Grating subunit stimulus

params.barWidth = 35;
iStim = ieStimulusGratingSubunit(params);
absorptions = iStim.absorptions;
movingBar = iStim;
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
os = osCreate('displayrgb');

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
x0 = -(numberElectrodesX-1)*electrodeArray.width/2;
y0 = -(numberElectrodesY-1)*electrodeArray.width/2;
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
        % electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/2)*(xPos-1) + electrodeArray.width, y0+(electrodeArray.width/2)*(yPos-1) + electrodeArray.width];
        electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/1)*(xPos-1) + 0, y0+(electrodeArray.width/1)*(yPos-1) + 0 + (mod(xPos,2)-.5)*(electrodeArray.width/2)];
    end
end
% xe = electrodeArray.center(:,:,1); ye = electrodeArray.center(:,:,2);
% figure; scatter(xe(:),ye(:));

th = (0:1/6:1)'*2*pi;
xh = electrodeArray.width/2*cos(th);
yh = electrodeArray.width/2*sin(th);
% figure
% fill(x,y,'r')
% axis square

% % Plot electrode array
eaSize = size(electrodeArray.center);
figure;
hold on;
for i = 1:eaSize(1)
    for j = 1:eaSize(2)
%         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
        plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
    end
end
axis equal
xlabel('Distance (m)'); ylabel('Distance (m)');
set(gca,'fontsize',14);

% Build the current stimulation activation window
% Gaussian activation from center of electrode
activationWindow = floor(retinalPatchSize(2)/numberElectrodesX);
electrodeArray.spatialWeight = fspecial('Gaussian', activationWindow, activationWindow/8);

% Visualize Gaussian activation
% figure; imagesc(electrodeArray.spatialWeight); 
% figure; surf(electrodeArray.spatialWeight); 
% xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
% title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);
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
            
%             sizeES = size(electrodeStimulus);
%             electrodeArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);

            % imagesc(electrodeStimulus); title(sprintf('%2.2f',mean(electrodeStimulus(:))));
        end
    end
end

% Show plot

% Apply Gaussian

%% Just electrode activation
% With < 20 uC/electrode there was no perception, Zrenner paper

% % Col activation
% electrodeArray.activation(:) = 0;
% electrodeArray.activation(4,:,4:8) = 0.1;
% electrodeArray.activation(4,:,24:28) = 0.5;
% electrodeArray.activation(4,:,44:48) = 0.75;
% electrodeArray.activation(4,:,64:68) = 1;

% name_str = 'col_electrode_10fps_lowOFF.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
%     for xPos = 2:2:numberElectrodesX-2
% %         for yPos = 2:2:numberElectrodesY
%             frame = frame + 20;
%             electrodeArray.activation(xPos,:,frame:frame+4) = 1;
%             
% %         end
%     end
%     electrodeArray.activation(xPos,yPos,frame+20) = 0;
% params.nSteps = frame+20;

%%% Single electrode activation
% name_str = 'single_electrode_10fps_lowOFF.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
%     for xPos = 2:2:numberElectrodesX
%         for yPos = 2:2:numberElectrodesY
%             frame = frame + 20;
%             electrodeArray.activation(xPos,yPos,frame:frame+4) = 1;
%             
%         end
%     end
%     electrodeArray.activation(xPos,yPos,frame+20) = 0;
% params.nSteps = frame+20;


% % % U activation
% name_str = 'U_10fps.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
% for rep = 1:3
%     electrodeArray.activation(3,1:4,frame+1:frame+4) = 1;
%     electrodeArray.activation(6,1:4,frame+1:frame+4) = 1;
%     electrodeArray.activation(3:6,4,frame+1:frame+4) = 1;
%     frame = frame+30;
% end


%%%%%%
% electrodeArray.activation(4,:,4:8) = 0.1;
% electrodeArray.activation(4,:,24:28) = 0.5;
% electrodeArray.activation(4,:,44:48) = 0.75;
% electrodeArray.activation(4,:,64:68) = 1;
%% Add 5 Hz spiking of stimulus

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 6;        % Radius in mm
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
hold on;
for spInd = 1:length(innerRetina.mosaic)
for i = 1:eaSize(1)
    for j = 1:eaSize(2)
        subplot(floor(sqrt(length(innerRetina.mosaic))),ceil(sqrt(length(innerRetina.mosaic))),spInd); 
        hold on;
%         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
        plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
    end
end
end

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

            end
        end
    end
end

innerRetinaInput = innerRetinaInput./10;%(numberElectrodesX*numberElectrodesY);

% figure; imagesc(squeeze(minXY(:,:,10,1,1)))
% figure; imagesc(squeeze(minXY(:,:,10,2,2)))
%% Build RGC activation functions

% figure; hold on;
% for i = 1:10
% % x0 = -.25;%*rand(1,1); 
% x0 = -.1 - .4*rand(1,1)
% thr = 10;%*rand(1,1);
% % figure; 
% xp = -.5:.01:1;
% plot((xp),0.5./(1+exp(-thr*(x0+(xp)))));
% end

% figure; hold on;
for mosaicInd = 1:length(innerRetina.mosaic)
    [yc xc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
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
    clear innerRetinaActivation  i0all xc yc
    [yc xc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
    i0all{mosaicInd} = .5*median(innerRetinaInput(:))*ones(xc,yc);% -.05 - .45*rand(xc,yc);
    for xind = 1:xc
        for yind = 1:yc
            
            for frame = 1:params.nSteps
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction{xind,yind,mosaicInd}(innerRetinaInput(xind,yind,mosaicInd));
                funcParams = innerRetinaThreshold{xind,yind,mosaicInd};
                thr = 80;%funcParams(1); 
                i0 = i0all{mosaicInd}(xind,yind);
                % innerRetinaFunction = @(iElectrode) 10./(1+exp(-thr*(i0+iElectrode)));
                % innerRetinaFunction = @(iElectrode) log(1./(1+exp(-thr*(i0+iElectrode))));
                
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction(innerRetinaInput(xind,yind,mosaicInd));
                
%                 innerRetinaFunction =  @(iElectrode) (50*(1./(1+exp(-thr*(i0+iElectrode)))));
%                 innerRetinaActivation{xind,yind}(frame) = innerRetinaFunction(innerRetinaInput(xind,yind,frame,mosaicInd));
                
                innerRetinaFunction = @(iElectrode) (5*iElectrode);
                innerRetinaActivation{xind,yind}(frame) = 50*innerRetinaFunction(innerRetinaInput(xind,yind,frame,mosaicInd));
            end
            mosaicSet(innerRetina.mosaic{mosaicInd},'responseLinear', innerRetinaActivation);
        end
    end
end
irPlot(innerRetina, 'linear');

irPlot(innerRetina,'mosaic');
% Visualize thresholds
% figure; hold on;
% for ii = 1:xc
%     for ji = 1:yc
%         xp = -.5:.01:1;
%         plot((xp),1./(1+exp(-thr*(i0all{mosaicInd}(ii,ji)+(xp)))));
%     end
% end
%% Compute RGC spiking
numberTrials = 1;
for tr = 1:numberTrials
    innerRetina = irComputeSpikes(innerRetina);
end
%% Invert representation to form image/movie
clear stimulusReconstruction
[stimulusReconstruction, paramsRec] = irReconstruct(innerRetina);


%% Build RGC array for healthy retina

clear paramsIR innerRetinaHealthy
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 5;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetinaHealthy = irCreate(os,paramsIR);
innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onMidget','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offMidget','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onParasol','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offParasol','model',model);

%%
innerRetinaHealthy = irComputeContinuous(innerRetinaHealthy,os);
numberTrials = 1;
for tr = 1:numberTrials
    innerRetinaHealthy = irComputeSpikes(innerRetinaHealthy);
end

irPlot(innerRetinaHealthy, 'linear');
% irPlot(innerRetinaHealthy, 'mosaic');
%% Invert representation to form image/movie
clear stimulusReconstructionHealthy
[stimulusReconstructionHealthy, paramsRecHealthy] = irReconstruct(innerRetinaHealthy);

%%
name_str = 'big_slow_grating2_20_withHealthy_OnMidget.mp4';
path_str = '/Users/james/Documents/MATLAB/isetbio misc/pixium_videos/';
vObj = VideoWriter([path_str name_str],'MPEG-4');
vObj.FrameRate = 10;
vObj.Quality = 100;
open(vObj);

% Play the movie with the stimulus
for loopv = 1%:10
h1=figure; set(gcf,'position',[160 60 1070 740]);
hold on;
for frame1 = 12:params.nSteps%size(movingBar.sceneRGB,3)
    subplot(221);
    imagesc(squeeze(movingBar.sceneRGB(:,:,frame1,:)));
    colormap gray; 
    
    subplot(222);
    for xPos = 1:numberElectrodesX
        for yPos = 1:numberElectrodesY
            hold on;
            fill(xh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,1),yh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,2),electrodeArray.activation(xPos,yPos,frame1))
        end
    end
    caxis([0 1]);
    
    
    subplot(223);    
    imagesc((stimulusReconstructionHealthy(1:paramsRecHealthy.maxx,1:paramsRecHealthy.maxy,frame1)));
     colormap gray
    caxis([.5*paramsRecHealthy.minR .5*paramsRecHealthy.maxR]);
%     caxis([0 .5*paramsRecHealthy.maxR]);
    title('Healthy');
    
    subplot(224);    
    imagesc(flipud(stimulusReconstruction(1:paramsRec.maxx,1:paramsRec.maxy,frame1)));
     colormap gray
%     caxis([.5*paramsRec.minR .5*paramsRec.maxR]);
    caxis([0 .5*paramsRec.maxR]);
    title('Prosthetic');
%     pause(0.1);
drawnow

    F = getframe(h1);
%     writeVideo(vObj,F);
end
end


close(vObj)