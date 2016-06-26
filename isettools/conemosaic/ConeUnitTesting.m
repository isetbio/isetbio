%% Conemosaic design and unit testing
%
%

%% Create scene
ieInit
scene = sceneCreate;

%% Create optical image
oi = oiCreate;
oi = oiCompute(oi,scene);

vcAddObject(oi); oiWindow;

%% Create cone mosaic
cMosaic = coneMosaic;
cMosaic.setSizeToFOV([1 0.8],'focalLength',oiGet(oi,'optics focallength'));


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
cMosaic.emPositions = cMosaic.emGenSequence(5000);
cMosaic.compute(oi);

%% Plot
cMosaic.plot('cone mosaic');
cMosaic.plot('cone fundamentals');
cMosaic.plot('macular transmittance');
cMosaic.plot('spectral qe');