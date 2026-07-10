% SkipFile
%% t_chessSet_cloud.m
%
% Render a very nice view of the chess set using the cloud. 
% 
% Depends on: iset3d, isetbio, Docker, isetcloud
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
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
saveDirName = sprintf('ChessSet_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic
dockerAccount= 'tlian';
dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://primal-surfer-140120.appspot.com';
clusterName = 'trisha-chess';
zone         = 'us-central1-a'; %'us-west1-a';    
instanceType = 'n1-highcpu-32';
gcp = gCloud('dockerAccount',dockerAccount,...
    'dockerImage',dockerImage,...
    'clusterName',clusterName,...
    'cloudBucket',cloudBucket,'zone',zone,'instanceType',instanceType);
toc
% gcp.Configlist; Doesn't seem to work

% Render depth?
gcp.renderDepth = true;

% Clear the target operations
gcp.targets = [];

%% Load the scene
myScene = sceneEye('chessSet');

%% Set parameters

myScene.accommodation = 1.43; 
myScene.pupilDiameter = 4;
myScene.fov = 30;

myScene.numCABands = 8;
myScene.diffractionEnabled = false;
myScene.numBounces = 4;

myScene.numRays = 1024;
myScene.resolution = 512;

myScene.eyePos = myScene.eyePos + [0.1 0.4 -0.3];
myScene.eyeTo = myScene.eyeTo + [0 0 -0.6];

% Zoom in 
forward = myScene.eyeTo - myScene.eyePos;
myScene.eyePos = myScene.eyePos + forward*0.2;

myScene.name = 'ChessSet';
    
% Normal render
% [oi, results] = myScene.render;
% ieAddObject(oi);
% oiWindow;

%% Render on cloud

uploadFlag = true;
[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,myScene,'uploadZip',uploadFlag);

%% Render
gcp.render();

%% Check for completion

% Save the gCloud object in case MATLAB closes
gCloudName = sprintf('%s_gcpBackup_%s',mfilename,currDate);
save(fullfile(myScene.workingDir,gCloudName),'gcp','saveDir');
    
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
    
    myScene = seAll{ii};
    
    saveFilename = fullfile(saveDir,[myScene.name '.mat']);
    save(saveFilename,'oi','myScene');
    
end

