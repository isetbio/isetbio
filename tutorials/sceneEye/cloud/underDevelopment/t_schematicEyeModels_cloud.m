% SkipFile
%% t_schematicEyeModels_cloud
% Render the same scene using a couple of different eye models.

%% Initialize
if piCamBio
    fprintf('%s: requires ISETBio, not ISETCam\n',mfilename); 
    return;
end
ieInit;
if ~mcDockerExists, mcDockerConfig; end % check whether we can use docker
if ~mcGcloudExists, mcGcloudConfig; end % check whether we can use google cloud sdk;

%% Initialize save folder
% Since rendering these images often takes a while, we will save out the
% optical images into a folder for later processing.
currDate = datestr(now,'mm-dd-yy_HH_MM');
saveDirName = sprintf('eyeModelComparison_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize cluster
tic

dockerAccount= 'tlian';
projectid = 'renderingfrl';
dockerImage = 'gcr.io/renderingfrl/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://renderingfrl';

clusterName = 'snellen';
zone         = 'us-central1-a';
instanceType = 'n1-highcpu-32';

gcp = gCloud('dockerAccount',dockerAccount,...
    'dockerImage',dockerImage,...
    'clusterName',clusterName,...
    'cloudBucket',cloudBucket,...
    'zone',zone,...
    'instanceType',instanceType,...
    'projectid',projectid,...
    'maxInstances',10);
toc

% Render depth
gcp.renderDepth = true;

% Clear the target operations
gcp.targets = [];

%% Load up a scene

thisScene = sceneEye('numbersAtDepth');

% Set general parameters
thisScene.fov = 30;
thisScene.resolution = 512;
thisScene.numRays = 2048;
thisScene.numCABands = 16;

%% Try the Navarro eye model

% This tell isetbio which model to use.
thisScene.modelName = 'Navarro';

% The Navarro model has accommodation, but let's set it to infinity for now
% since other models may not have accommodation modeling.
thisScene.accommodation = 0;

% Upload to cloud
thisScene.name = 'navarro'; % The name of the optical image
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,thisScene,'uploadZip',false);

%% Try the Gullstrand-LeGrand Model

% The gullstrand has no accommodation modeling. 
thisScene.modelName = 'Gullstrand';

% Upload to cloud
thisScene.name = 'gullstrand'; % The name of the optical image
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,thisScene,'uploadZip',false);

%% Try Arizona eye model

thisScene.modelName = 'Arizona';
thisScene.accommodation = 0;

% Upload to cloud
thisScene.name = 'arizona'; % The name of the optical image
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,thisScene,'uploadZip',true);

%% Compare Navarro and Arizona eye accommodated

thisScene.modelName = 'Navarro';
thisScene.accommodation = 5;

% Upload to cloud
thisScene.name = 'navarroAccommodated'; % The name of the optical image
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,thisScene,'uploadZip',false);

thisScene.modelName = 'Arizona';
thisScene.accommodation = 5;

% Upload to cloud
thisScene.name = 'arizonaAccommodated'; % The name of the optical image
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,thisScene,'uploadZip',true);

%% Render
gcp.render();

%% Check for completion

% Save the gCloud object in case MATLAB closes
gCloudName = sprintf('%s_gcpBackup_%s',mfilename,currDate);
save(fullfile(saveDir,gCloudName),'gcp','saveDir');

% Pause for user input (wait until gCloud job is done)
x = 'N';
while(~strcmp(x,'Y'))
    x = input('Did the gCloud render finish yet? (Y/N)','s');
end


%% Download the data

[oiAll, seAll] = downloadFromCloud(gcp);

for ii=1:length(oiAll)
    
    oi = oiAll{ii};
    
    ieAddObject(oi);
    oiWindow;
    
    thisScene = seAll{ii};
    
    saveFilename = fullfile(saveDir,[thisScene.name '.mat']);
    save(saveFilename,'oi','thisScene');
    
end
