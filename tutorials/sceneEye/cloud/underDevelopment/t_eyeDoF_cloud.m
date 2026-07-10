% SkipFile
%% t_eyeDoF.m
%
% This tutorial shows the effect of pupil diameter on the depth of field in
% the scene.
% 
% Depends on: iset3d, isetbio, Docker, isetcloud
%
% TL ISETBIO Team, 2017 

%% Initialize ISETBIO
if isequal(piCamBio,'isetcam')
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
saveDirName = sprintf('dof_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic
dockerAccount= 'tlian';
dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://primal-surfer-140120.appspot.com';
clusterName = 'tl-dof';
zone         = 'us-central1-b'; %'us-west1-a';    
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

% The "chessSetScaled" is the chessSet scene but scaled and shifted in a
% way that emphasizes the depth of field of the eye. The size of the chess
% pieces and the board may no longer match the real world.
myScene = sceneEye('chessSetScaled');

%% Set fixed parameters
myScene.accommodation = 1/0.28; 
myScene.fov = 30;

myScene.numCABands = 6;
myScene.diffractionEnabled = false;
myScene.numBounces = 3;

lqFlag = false;

%% Loop and upload to cloud

pupilDiameters = [2 2.5 3 3.5 4 4.5 5 5.5 6];
% pupilDiameters = [2 4 6];

for ii = 1:length(pupilDiameters)
    
    myScene.pupilDiameter = pupilDiameters(ii);
    
    % Change number of rays depending on pupil size
    if(~lqFlag)
        if(pupilDiameters(ii) < 3)
            myScene.numRays = 8192;
        elseif(pupilDiameters(ii) >= 4)
            myScene.numRays = 2048;
        else
            myScene.numRays = 4096;
        end
    else
        myScene.numRays = 256;
    end
    
    if(lqFlag)
        myScene.resolution = 256;
    else
        myScene.resolution = 512;
    end
    
    myScene.name = sprintf('DoF%0.2fmm',pupilDiameters(ii));

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





