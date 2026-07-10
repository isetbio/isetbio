% SkipFile
%% t_accommodationMTF_fineTune.m
%
% Render a slanted bar through space. Accommodate always at the plane.
% Measure the MTF as the plane moves. We would expect the curves to always
% be the same. 
%
% This can help verify the accommodation modeling. 
%
% Fine tune - see if the accommodation is just slightly off. If so, how off
% is it?
%
% Depends on: iset3d, isetbio, Docker, isetcloud
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;
if ~mcDockerExists, mcDockerConfig; end % check whether we can use docker
if ~mcGcloudExists, mcGcloudConfig; end % check whether we can use google cloud sdk;

%% Initialize save folder
% Since rendering these images often takes a while, we will save out the
% optical images into a folder for later processing.
currDate = datestr(now,'mm-dd-yy_HH_MM');
saveDirName = sprintf('accomVerFT_%s',currDate);
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

clusterName = 'accom-ver';
zone         = 'us-west1-a';    
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

%% Loop through accommodations

accommodations = 6.5; %linspace(9.62,9.75,5);
distToPlane = 1/6.5
    
% Load scene with plane at a specific distance
myScene = sceneEye('slantedBar','planeDistance',distToPlane);

nRenders = length(accommodations);

for ii = 1:nRenders
    
    currAccomm = accommodations(ii);
 
    myScene.name = sprintf('mtf_accom_%0.2fdp',currAccomm);
    
    myScene.accommodation = currAccomm;
    myScene.pupilDiameter = 3;
    
    myScene.fov = 1;
    
    % HQ version
%     myScene.numCABands = 16;
%     myScene.numRays = 512;
%     myScene.resolution = 512;

    % LQ version
    myScene.numCABands = 0;
    myScene.numRays = 256;
    myScene.resolution = 128;
    
%     if(ii == nRenders)
%         fprintf('Uploading zip... \n');
%         uploadFlag = true;
%     else
%         uploadFlag = false;
%     end
%     [cloudFolder,zipFileName] =  ...
%         sendToCloud(gcp,myScene,'uploadZip',uploadFlag);

    % Normal render
    [oi, results] = myScene.render();
    ieAddObject(oi);
    oiWindow;
    
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

