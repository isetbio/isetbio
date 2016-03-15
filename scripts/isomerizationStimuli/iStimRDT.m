%% Script shows how we store iStim up in the RDT
%
%




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

%% Create the file with the bar iStim

clear params
params.barWidth = 10; 
params.fov=0.4;
params.timeInterval = 0.002;  % Two ms time interval
iStim = ieStimulusBar(params);

fname = fullfile(isetbioRootPath,'local','barMovie.mat');
save(fname,'iStim');

%%
rd.publishArtifact(fname);

%%