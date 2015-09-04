%% t_VernierSimonsSlides
%
% This tutorial script uses biological and computational methods to analyze
% vernier acuity (super acuity) in human vision
%
% Vernier acuity (or positional acuity) is a measurement of sensitivity of
% human eye in detecting mis-alignment of simple object (lines, etc.)
%
% In this script, we compute the irradiance and human cone absorptions for
% a scene with two mis-aligned lines. Then, we try to discriminate the
% aligned and mis-aligned cases by using first order statistics and machine
% learning classifiers
%
%  HJ/BW, ISETBIO TEAM, 2015

%  Initialize a new session
ieInit;

%% Create the display

% In this example we impose a linear gamma table, though
% in general it could be the default or anything.
dpi = 200; d = displayCreate('LCD-Apple','dpi',dpi);
% d = displaySet(d, 'gamma', repmat(linspace(0, 1, 256)', [1 3]));

viewDist = 2; % viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);
d = displaySet(d,'g table','linear');

%% Create the the RGB image
[~, p] = imageVernier();   % Mainly to get the parameters
p.pattern = 0*ones(1,65); p.pattern(33) = 1;

% Aligned
p.offset = 0; imgA = imageVernier(p);
% vcNewGraphWin; imagesc(imgA)

% Misaligned
p.offset = 2; imgM = imageVernier(p);
% vcNewGraphWin; imagesc(imgM)

meanLum = 50; 
sceneA = sceneFromFile(imgA, 'rgb', meanLum, d); % aligned
sceneM = sceneFromFile(imgM, 'rgb', meanLum, d); % mis-aligned

vcAddObject(sceneA); 
vcAddObject(sceneM); 
sceneWindow;
%% Examine some of the scene properties

% This is the scene luminance at the different sample points on the display
sz = sceneGet(sceneM,'size');
scenePlot(sceneM,'luminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
scenePlot(sceneM,'radiance hline',[1,round(sz(1)/2)]);

%% Compute Irradiance with Optics Wavefront
%    In this section, we compute the optical image (irradiance map) by
%    using human optics wavefront.
%
%    If we just want a standard human optics model, we can simplify the
%    code as oiCraete('wvf human');

%  Load Zernike coefficient
pupilSize = 3; % pupil size in mm
zCoefs = wvfLoadThibosVirtualEyes(pupilSize);

%  Create wavefront structure
wave= sceneGet(sceneA,'wave');
wvf = wvfCreate('wave', wave, 'zcoeffs', zCoefs, 'name', 'human optics');
wvf = wvfSet(wvf, 'calc pupil size', pupilSize); 

% Adjust for individuals
% Here, we use defocus as an example. For more adjustable entries, see
% wvfOSAIndexToVectorIndex

% ajust zernike coefficient for defocus 
% if we need to use defocus in diopters, use wvfDefocusDioptersToMicrons to
% do the conversion
zDefocus = -0.0104; 
wvf = wvfSet(wvf, 'zcoeffs', zDefocus, {'defocus'});

% compute psf and convert to optical image structure
wvf = wvfComputePSF(wvf);
oi = wvf2oi(wvf, 'human');
oi = oiSet(oi, 'name', sprintf('Human WVF %.1f mm', pupilSize));

%% Examine the point spread function at different wavelengths

% This is the standard model based on the Thibos AO estimates of the human
% wavefront function
oiPlot(oi,'psf 550');
nPoints = 50;
oiPlot(oi,'ls wavelength',[],nPoints);

%% compute irradiance map (optical image)
oiA = oiCompute(sceneA, oi);
oiM = oiCompute(sceneM, oi);
rect = [9     8    64    65];
oiA = oiCrop(oiA,rect);
oiM = oiCrop(oiM,rect);
vcAddObject(oiA); vcAddObject(oiM); oiWindow;

% Another way to do this computation is using the chromatic aberration in
% the Marimont and Wandell model (1994, JOSA).  You can create that oi
% structure simply by calling oi = oiCreate('human');

%% Examine the irradiance at the retina, prior to absorption

% This is the scene luminance at the different sample points on the display
sz = sceneGet(oiA,'size');
oiPlot(oiA,'illuminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
oiPlot(oiA,'irradiance hline',[1,round(sz(1)/2)]);
view(135,23);
% Notice that the short-wavelength light is spread a lot more at the
% retinal surface than the middle wavelength light.
% The long-wavelength is spread a little bit more.

%% Compute Photon Absorptions of Human Cone Receptors

% Compute the human cone absorption samples with fixational eye movement.
cones = sensorCreate('human');
cones = sensorSetSizeToFOV(cones, sceneGet(sceneA, 'fov'), sceneA, oiA);
cones = sensorCompute(cones,oiM);
vcAddObject(cones); sensorWindow('scale',1);

%% Static analysis, without eye movements

% Plot the average of the top and bottom half, looking for displacement
% of peak
e  = sensorGet(cones,'photons');
sz = sensorGet(cones,'size');
midSensor = floor(sz(2)/2) + 1;

midRow    = floor((sz(1)/2));
topSensor = 1:midRow;
botSensor = (midRow+1):sz(1);
topE  = mean(e(topSensor,:),1)/max(e(:));
botE  = mean(e(botSensor,:),1)/max(e(:));

% Graph the top and bottom
vcNewGraphWin;
plot([topE(:),botE(:)])
% line([midSensor,midSensor],[min(e(:)) max(e(:))],'color','r')
grid on; legend('top','bottom')
title(sprintf('Bar offset %0.f (arc sec)',sceneGet(sceneM,'degrees per sample','arcsec')))
ylabel('Normalized cone absorptions');
xlabel('Position (um)');

%% Analysis with eye movements
nFrames = 50;        % Number of exposures samples

expTime = sensorGet(cones, 'exp time');   % Usually about 50 ms
emDuration = 0.001;
emPerExposure = expTime / emDuration;
sensor = sensorSet(cones, 'exp time', emDuration);

% Generate eyemovement
p.nSamples = nFrames * 50;
p.emFlag   = [1 1 1];       % Include tremor drift and saccade
cones = eyemoveInit(cones, p);

% % We could enlarge the eye movement tremor
% tremor = sensorGet(cones,'em tremor');
% tremor.amplitude = tremor.amplitude*2;
% cones = sensorSet(cones,'em tremor',tremor);

% Generate the eye movement sequence again
% cones = emGenSequence(cones);

% The coneAbsorptions function is an interface to sensorCompute. Notice
% that when we make an eye movement video we call coneAbsorptions, not
% sensorCompute.
cones = coneAbsorptions(cones, oiM);  % Use oiA for aligned

% Show the eye movement positions
ePos = sensorGet(cones,'sensor positions');
vcNewGraphWin;
plot(ePos(:,1),ePos(:,2),'-o')
xlabel('Cone position');
ylabel('Cone Position')
% set(gca,'xlim',[-15 15],'ylim',[-15 15]);
grid on

fprintf('%d eye positions\n',size(ePos,1))
fprintf('%d integration periods\n',nFrames);

%% Show the cone absorptions each millisecond

% Get the time series out from the cone photon data
absorptions = sensorGet(cones,'photons');
absorptions = ieScale(absorptions,150);

% This should work, but it depends on having a later version of Matlab than
% 2013b
%   mplay(absorptions,'intensity',10);

%% Movie of the cone absorptions

step = 20;
tmp = coneImageActivity(cones,[],step,false);

% Show the movie
vcNewGraphWin;
tmp = tmp/max(tmp(:));
for ii=1:size(tmp,4)
    img = squeeze(tmp(:,:,:,ii));
    imshow(img.^0.3); truesize;
    title('Cone absorptions')
    drawnow
end

%% Black and white version
vcNewGraphWin;
% vObj = VideoWriter('coneAbsorptions.avi');
% open(vObj);
colormap(gray);
nframes = size(absorptions,3);
% Record the movie
for j = 1:step:nframes 
    image(absorptions(:,:,j)); drawnow;
    title('Cone absorptions')
%     F = getframe;
%     writegit puVideo(vObj,F);
end
% close(vObj);
fprintf('Max cone absorptions %.0f\n',max(absorptions(:)));

%% Temporal dynamics applied to the cone absorptions

% [cones,adaptedData] = coneAdaptAlt(cones,'rieke');
[adaptedData, adaptedOS] = coneAdaptAlt(cones,'rieke');


% The adapted data values are all negative current.
adaptedData = ieScale(adaptedData,0,1);

step = 20;
tmp = coneImageActivity(cones,adaptedData,step,false);

% Show the movie
vcNewGraphWin;
tmp = tmp/max(tmp(:));
for ii=1:size(tmp,4)
    img = squeeze(tmp(:,:,:,ii));
    imshow(img.^0.3); truesize;
    title('Cone photocurrent')
    drawnow
end


%% Black and white version
vcNewGraphWin;
% vObj = VideoWriter('coneVoltage.avi');
%  open(vObj);
adaptedData = 150*ieScale(adaptedData,0,1);
colormap(gray);
nframes = size(adaptedData,3);
% Record the movie
for j = 1:step:nframes
    image(adaptedData(:,:,j));
    title('Cone photocurrent');
    drawnow;
    %     F = getframe;
    %     writeVideo(vObj,F);
end
%  close(vObj);



%% RGC Responses

%% Create RGC structure

% creating the RGC parameter object
rgcP = rgcParameters;

% The absorptions structure is also the RGC parameter data
% rgcP.set('scene',scene);
rgcP.set('oi',oi);
rgcP.set('sensor',cones);
rgcP.set('cone voltages',adaptedOS.ConeCurrentSignal); % check units on cone voltages

% What is the default?  RGC spacing the same as cone spacing (I think).
rgcP.addLayer('on parasol', 20);  
rgcP.addLayer('off parasol');  
 
rgcP.addLayer('on midget');  
rgcP.addLayer('off midget');  

rgcP.addLayer('small bistratified');  % Sign of signal seems wrong. Check.

%%
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
% figure; imagesc(rgcP.layers{1,1}.currentSpkTS')


%%
figure; 
for spind = 1:5
    subplot(3,2,spind); 
%     imagesc(rgcP.layers{1,spind}.currentSpkTS);
        imagesc(rgcP.layers{1,spind}.rgcvTimeSeries);
    switch spind
        case 1
            title('on parasol');
        case 2
            title('off parsol');
        case 3
            title('on midget');
        case 4 
            title('off midget');
        case 5
            title('small bistratified');
    end
xlabel('time (ms)'); ylabel('RGC Spikes');
end

%% Black and white version
vcNewGraphWin;
% vObj = VideoWriter('coneVoltage.avi');
%  open(vObj);
% adaptedData = 150*ieScale(adaptedData,0,1);
spind = 2;
rsv = rgcP.layers{1,spind}.gridSize;
adaptedData = 150*ieScale(reshape(rgcP.layers{1,spind}.rgcvTimeSeries,rsv(1),rsv(2),2500),0,1);
% adaptedData = 150*ieScale(reshape(rgcP.layers{1,spind}.currentSpkTS,rsv(1),rsv(2),2500),0,1);
colormap('gray');
nframes = size(adaptedData,3);
% Record the movie
for j = 1:step:nframes
    image(adaptedData(:,:,j));
    switch spind
        case 1
            title('on parasol');
        case 2
            title('off parsol');
        case 3
            title('on midget');
        case 4 
            title('off midget');
        case 5
            title('small bistratified');
    end
    drawnow;
    %     F = getframe;
    %     writeVideo(vObj,F);
end
%  close(vObj