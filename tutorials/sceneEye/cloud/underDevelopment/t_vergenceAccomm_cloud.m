% SkipFile
%% t_vergenceAccomm.m
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
saveDirName = sprintf('va_%s',currDate);
saveDir = fullfile(isetbioRootPath,'local',saveDirName);
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% Initialize your cluster
tic
dockerAccount= 'tlian';
dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v3-spectral-gcloud';
cloudBucket = 'gs://primal-surfer-140120.appspot.com';
clusterName = 'trisha-va';
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
myScene = sceneEye('chessSet');

%% Set parameters

% LQ
myScene.resolution = 128; 
myScene.numRays = 256;
myScene.numCABands = 0;

% HQ
% myScene.resolution = 512; 
% myScene.numRays = 2048;
% myScene.numCABands = 6;

ipd = 64e-3; % Average interpupillary distance

originalPos = myScene.eyePos;
originalWorld = myScene.recipe.world;
startingPos = originalPos + [0 0 0.1];

%% Create binocular retinal images

% vergenceZ = [-0.2 0.0 0.2]; % Along z-axis
vergenceZ = linspace(-0.2,0.2,10);

for ii = 1:length(vergenceZ)
    
    vergencePoint = [0 startingPos(2) vergenceZ(ii)];
    
    leftEyePos = startingPos - [ipd/2 0 0];
    rightEyePos = startingPos + [ipd/2 0 0];
    
    myScene.eyeTo = vergencePoint;
    
    % Reset world (Important!)
    myScene.recipe.world = originalWorld;
    
    % Add the target sphere
    myScene.recipe = piAddSphere(myScene.recipe,...
        'rgb',[1 0 0],...
        'radius',0.005,...
        'location',vergencePoint);
    
    % Set accommodation to the right distance
    dist = sqrt(sum((vergencePoint -leftEyePos).^2)); % in mm
    myScene.accommodation = 1/dist;
    
    % Gcloud doesn't like the "-" sign
    vergenceStr = num2str(vergenceZ(ii));
    vergenceStr = strrep(vergenceStr,'-','neg');
    
    % Left eye
    myScene.eyePos = leftEyePos;
    myScene.name = sprintf('leftEye_%sm',vergenceStr);
    sendToCloud(gcp,myScene,'uploadZip',false); 
    
    % Right Eye
    myScene.eyePos = rightEyePos;
    myScene.name = sprintf('rightEye_%sm',vergenceStr);
    % Upload zip for final image
    if(ii == length(vergenceZ))
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
