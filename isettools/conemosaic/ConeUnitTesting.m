%% Conemosaic design and unit testing
%
%

%% Create scene
ieInit
scene = sceneCreate('slanted bar');

%% Create optical image
oi = oiCreate;
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;

% oi = oiSet(oi, 'lens density', 0);
% oi = oiCompute(oi,scene);
% vcAddObject(oi); oiWindow;

%% Create cone mosaic
cMosaic = coneMosaic('os', osBioPhys);

%% Comparing with sensor calculation
sensor = sensorCreate;
sensor = sensorSet(sensor,'pattern',cMosaic.pattern);
sensor = sensorSet(sensor,'noise flag',0);
sensor = sensorCompute(sensor,oi);
photons = sensorGet(sensor,'photons');
vcNewGraphWin; imagesc(photons);

cMosaic.noiseFlag = 0;
cMosaic.compute(oi);
vcNewGraphWin; imagesc(cMosaic.absorptions);

vcNewGraphWin; plot(photons(:), cMosaic.absorptions(:),'.');
identityLine;

vcNewGraphWin;
hist(photons(:) -  cMosaic.absorptions(:),100);

%% Test noise addition
cMosaic.noiseFlag = true;
cMosaic.compute(oi);
nAbsorptions = cMosaic.absorptions;
vcNewGraphWin; imagesc(nAbsorptions);

cMosaic.noiseFlag = false;
cMosaic.compute(oi);
mAbsorptions = cMosaic.absorptions;
vcNewGraphWin; hist(nAbsorptions(:)-mAbsorptions(:), 100);

%% Eye movement testing
cMosaic.emGenSequence(5000);
cMosaic.compute(oi);
m = mean(cMosaic.absorptions,3);
vcNewGraphWin; imagesc(m)

vcNewGraphWin;
plot(cMosaic.emPositions(:,1),cMosaic.emPositions(:,2),'-');

%% Plot
cMosaic.plot('cone mosaic');
cMosaic.plot('cone fundamentals');
cMosaic.plot('macular transmittance');
cMosaic.plot('eye spectral qe','oi',oi);

%% Show the time series of absorptions as a movie
[~,uData] = cMosaic.plot('absorptions');

%% Compute photocurrent
[~, current] = cMosaic.compute(oi, 'currentFlag', true);

%% GUI
cMosaic.guiWindow;