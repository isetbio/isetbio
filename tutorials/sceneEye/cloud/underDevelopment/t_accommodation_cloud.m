% SkipFile
%% t_accommodation_cloud.m
%
% Demonstrate sceneEye accommodation using the numbersAtDepth scene and
% rendering using the cloud.

% Depends on: iset3d, isetbio, Docker, isetcloud.
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
saveDirName = sprintf('accom_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic
dockerAccount= 'tlian';
dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://primal-surfer-140120.appspot.com';
clusterName = 'trisha-accom';
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

%% Select scene

myScene = sceneEye('numbersAtDepth');

myScene.eyePos = myScene.eyePos + [0 0.005 0];
myScene.fov = 35;

%% Set fixed parameters

myScene.numRays = 256;
myScene.resolution = 256;

myScene.pupilDiameter = 4;
myScene.numCABands = 6;

%% Step through accommodation

accomm = [3:10]; % in diopters

for ii = 1:length(accomm)
    
    myScene.accommodation = accomm(ii);
    myScene.name = sprintf('accom_%0.2fdpt',myScene.accommodation);
    
    % Instead of rendering, we add it to the list of gcloud targets
    % Note: since sceneEye is a a special object (different from the usual
    % recipe we use in iset3d) we use a helper function that will
    % facilitate the adding of each sceneEye into the gCloud class
    if(ii == length(accomm))
        % Upload the zip file during the final loop.  
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



