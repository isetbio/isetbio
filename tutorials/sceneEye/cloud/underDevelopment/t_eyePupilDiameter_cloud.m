% SkipFile
%% t_eyePupilDiameter_cloud.m
%
% Render the slanted bar at different pupil diameters.
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
saveDirName = sprintf('pupil_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic
dockerAccount= 'tlian';
dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://primal-surfer-140120.appspot.com';
clusterName = 'trisha-pupil';
zone         = 'us-central1-a'; %'us-west1-a';    
instanceType = 'n1-highcpu-32';
gcp = gCloud('dockerAccount',dockerAccount,...
    'dockerImage',dockerImage,...
    'clusterName',clusterName,...
    'cloudBucket',cloudBucket,'zone',zone,'instanceType',instanceType);
toc

% Render depth
gcp.renderDepth = true;

% Clear the target operations
gcp.targets = [];

%% Load scene
planeDistance = 5; % 5 meters away
myScene = sceneEye('slantedBar','planeDistance',planeDistance); 

%% Set fixed parameters
myScene.accommodation = 1/planeDistance; % Accomodate to plane
myScene.fov = 2;

myScene.numCABands = 8;
myScene.diffractionEnabled = true;

%% Loop through pupil diameters

% HQ version
% myScene.resolution = 1024; 
% pupilDiameters = [6 5 4 3 2 1];
% numRays = [2048 2048 2048 4096 8192 8192]; % This needs to change with pupil diameter

% Fast test version
pupilDiameters = [6 5 4 3 2];
numRays = [128 128 128 128 128];
myScene.resolution = 128; 

if(length(numRays) ~= length(pupilDiameters))
    error('numRays and pupilDiameters length need to match!')
end

for ii = 1:length(pupilDiameters)
    
    currPupilDiam = pupilDiameters(ii);
    currNumRays = numRays(ii);
    
    myScene.pupilDiameter = currPupilDiam;
    myScene.numRays = currNumRays;
    
    myScene.name = sprintf('pupilDiam_%0.2fmm',currPupilDiam);
    
    if(ii == length(pupilDiameters))
        fprintf('Uploading zip... \n');
        uploadFlag = true;
    else
        uploadFlag = false;
    end
    [cloudFolder,zipFileName] =  ...
        sendToCloud(gcp,myScene,'uploadZip',uploadFlag);
    
end

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
    
    myScene = seAll{ii};
    
    saveFilename = fullfile(saveDir,[myScene.name '.mat']);
    save(saveFilename,'oi','myScene');
    
end




