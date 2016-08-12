%% Store iStim data in the RDT
%
% These files are used for testing RGC ... other purposes?
%
% JRG/BW ISETBIO Team, 2016


%% Create the moving bar iStim for RDT upload

clear params
params.barWidth = 10; 
params.fov=0.4;
params.timeInterval = 0.002;  % Two ms time interval
iStim = ieStimulusBar(params);

fname = fullfile(isetbioRootPath,'local','barMovie.mat');
save(fname,'iStim');


%% Create the file with the gabor movie iStim in it
params.freq = 6;
params.nSteps = 50;
params.GaborFlag = 0.2;
iStim = ieStimulusGabor(params); %#ok<NASGU>

fname = fullfile(isetbioRootPath,'local','gaborMovie.mat');
save(fname,'iStim');

%% Open the RDT
rd = RdtClient('isetbio');
rd.credentialsDialog;
rd.crp('/resources/data/istim');

%%
rd.publishArtifact(fname);


%%
rd.publishArtifact(fname);

%%