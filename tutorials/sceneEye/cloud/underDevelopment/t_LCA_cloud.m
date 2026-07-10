% SkipFile
%% t_LCA_cloud.m
%
% Render a slanted bar while moving the retina along the optical axis. We
% can use the images generated here to calculate the amount of LCA present
% in the eye
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
saveDirName = sprintf('LCA_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic

dockerAccount= 'tlian';
projectid = 'renderingfrl';
dockerImage = 'gcr.io/renderingfrl/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://renderingfrl';

clusterName = 'lca';
zone         = 'us-central1-a';    
instanceType = 'n1-highcpu-32';

gcp = gCloud('dockerAccount',dockerAccount,...
    'dockerImage',dockerImage,...
    'clusterName',clusterName,...
    'cloudBucket',cloudBucket,...
    'zone',zone,...
    'instanceType',instanceType,...
    'projectid',projectid);
toc

% Render depth
gcp.renderDepth = true;

% Clear the target operations
gcp.targets = [];

%% Turn on chromatic aberration to show color fringing.

distToPlane = 0.2;

% Move the retina plane 
retinaDistance = 16.00:0.05:16.60;
% retinaDistance = [16.3 16.2];
retinaRadius = 10000; % Make it flat
retinaSemiDiam = 0.15;

for ii = 1:length(retinaDistance)
    
    % Load scene with plane at a specific distance
    myScene = sceneEye('slantedBar','planeDistance',distToPlane);
    
    myScene.name = sprintf('slantedBar_LCA_%0.3fmm',retinaDistance(ii));
    
    % Calculate FOV needed to get the same image size
    myScene.fov = 2*atand(retinaSemiDiam/retinaDistance(ii));
    
    myScene.retinaDistance = retinaDistance(ii);
    myScene.accommodation = 1/distToPlane;
    myScene.pupilDiameter = 4;

    myScene.numCABands = 16;
    myScene.numRays = 1024;
    myScene.resolution = 512;
%    myScene.numRays = 128;
%    myScene.resolution = 128;
    
    if(ii == length(retinaDistance))
        fprintf('Uploading zip... \n');
        uploadFlag = true;
    else
        uploadFlag = false;
    end
    [cloudFolder,zipFileName] =  ...
        sendToCloud(gcp,myScene,'uploadZip',uploadFlag);

    % Normal render
%     [oi, results] = myScene.render();
%     ieAddObject(oi);
%     oiWindow;
    
end

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

