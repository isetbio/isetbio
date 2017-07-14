%%  Make an interesting cone mosaic data set and put it up on RDT
%
% s_cmRDT
%
% BW ISETBIO Team, 2017

%% Vernier example
% Needs WL/WLVernierAcuity
%

% Show the dependence on the cone mosaic size for the computational
% observer.
nTrials = 10;
nBasis  = 40;

% Integration time
tStep   = 10;         % Adequate for photocurrent (ms)

% Cone mosaic field of view in degrees
coneMosaicFOV = 0.35;

% Original scene
sceneFOV = 0.4;

% Spatial scale to control visual angle of each display pixel The rule is
% 6/sc arc sec for a 0.35 deg scene. If you change the scene to 0.5 deg
% then 0.5/0.35
sc = 3*(sceneFOV/0.35);   % If you do not multiply by a scalar, offset is 6 arc sec

s_EIParameters;
params.em.emFlag = [0 0 0]';

% Make the bar length a little less than the scene size
params.vernier.barLength = params.vernier.sceneSz(1)-1;
params.tsamples  = (-200:tStep:400)*1e-3;
%% Read the stimulus, which might have been saved

params.vernier.offset = 0;
[aligned, offset, ~, ~] = vaStimuli(params);

%%  Compute absorptions for multiple trials
tSamples = aligned.length;
cMosaic = coneMosaic;

% Sometimes we set the mosaic size to 15 minutes (.25 deg) because that is
% the spatial pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(params.cmFOV);

% Not sure why these have to match, but there is a bug if they don't.
cMosaic.integrationTime = aligned.timeStep;
cMosaic.noiseFlag = 'random';
% cMosaic.plot('impulse response'); cMosaic.plot('os impulse response');

cMosaic.pattern = 4*ones(size(cMosaic.pattern));
%% For aligned or offset

% Turn off eye movements for a moment


disp('Computing cone mosaic eye movements');
emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
    'em', params.em);

% compute absorptions for aligned and offset
disp('Computing cone mosaic current');
[~,alignedC] = cMosaic.compute(aligned, 'currentFlag', true, ...
    'emPaths', emPaths);

% Have a look cMosaic.window;

coneMosaicFile = fullfile(isetbioRootPath,'local','coneMosaicDataScones.mat');
save(coneMosaicFile, 'cMosaic','alignedC');

%%

rdt = RdtClient('isetbio');
rdt.credentialsDialog();  % wandell, Jxxx4XX

% You can change other places.
rdt.crp('/resources/data/cmosaics')
version1 = '1';
rdt.publishArtifact(coneMosaicFile, 'version', version1);

%% Make a cone mosaic of a face image

scene = sceneFromFile;
scene = sceneSet(scene,'fov',2);
scene = sceneAdjustIlluminant(scene,'D65');
vcAddObject(scene); sceneWindow;

oi = oiCreate;
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(sceneGet(scene,'fov'));
emPaths  = cMosaic.emGenSequence(60);

% Better way to compute, like above.  Does it require multiple trials?
cMosaic.compute(oi, 'emPaths', emPaths);
cMosaic.computeCurrent;
alignedC = cMosaic.current;

cMosaic.window;
faceFile = fullfile(isetbioRootPath,'local','coneMosaicDataF.mat');
save(faceFile,'cMosaic','alignedC');

rdt = RdtClient('isetbio');
rdt.credentialsDialog();  % wandell, Jxxx4XX

% You can change other places.
rdt.crp('/resources/data/cmosaics')
rdt.listArtifacts('recurseive',true,'print',true);
version1 = '1';
rdt.publishArtifact(faceFile, 'version', version1);

%%

