%% Conemosaic design and unit testing
%
%

%%
ieInit
scene = sceneCreate;

%% The new oi has the human lens in it.

oi = oiCreate;
oi = oiCompute(oi,scene);

vcAddObject(oi); oiWindow;

% TODO:  Make lensGet/lensSet for the new lens class
%%
cMosaic     = coneMosaic;

%%
% absorptions = cMosaic.computeSingleFrame(oi);

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

vcNewGraphWin; plot(photons(:),absorptions(:),'.');
identityLine;

%%  This routine adds noise

cMosaic.noiseFlag = true;
cMosaic.compute(oi);
nAbsorptions = cMosaic.absorptions;
vcNewGraphWin; imagesc(nAbsorptions);

cMosaic.noiseFlag = false;
cMosaic.compute(oi);
mAbsorptions = cMosaic.absorptions;
vcNewGraphWin; hist(nAbsorptions(:)-mAbsorptions(:),100);

%% Eye movement testing

cMosaic.emPositions = round(10*rand(10,2));

