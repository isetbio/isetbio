% SkipFile
%% t_PSFoverDefocus.m
%
% Depends on: iset3d, isetbio, Docker, RemoteDataToolbox
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
saveDirName = sprintf('defocusPSF_%s',currDate);
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

clusterName = 'defocus';
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


%% Generate a PSF at focus using the eye model

% Load scene
% A black disk located 100 meters away (on the y-axis) and spanning 120 deg FOV.
myScene = sceneEye('blackBackdrop');

pixelFOV = 0.0020; % Essentially a point. See calculation in s_comparePSF.m

pointDistance = 1/5; % In meters
pointWidth = 2*pointDistance*tand(pixelFOV/2);

% Add equal energy sphere on the optical axis. We will make it the same size as the
% pixel above.
myScene.recipe = piAddSphere(myScene.recipe,...
    'spectrum',[400 1 800 1],...
    'radius',pointWidth/2,...
    'location',[0 pointDistance 0]);

% Set rendering parameters
% Drop FOV to have higher chance of hitting the point
fovScale = 0.2;
oiFOV = 1.2500; % From s_comparePSF.m
myScene.fov = oiFOV*fovScale;

% Calculate new resolution
oiSS = 5.7952e-07; % From s_comparePSF.m
myScene.resolution = round(myScene.width/(oiSS*10^3));
myScene.numRays = 2^15; % This needs to be high, since the point is small and hard to hit!  
% scene3d.numRays = 256; % Test for now

% Add chromatic aberration
myScene.numCABands = 16;

% Set pupil diameter
myScene.pupilDiameter = 4;

pedestalDefocus = 1/pointDistance;
deltaDefocus = [-1 -0.5 0 0.5 1];
% deltaDefocus = [-1 0 1]; % Test

nRenders = length(deltaDefocus);

% Loop over accommodations
for ii = 1:nRenders
    
    % Change accommodation
    myScene.accommodation = pedestalDefocus + deltaDefocus(ii);
    
    % Set name
    myScene.name = sprintf('psf_3deye_%0.2fdpt',myScene.accommodation);
    
    % Render with cloud
    if(ii == nRenders)
        fprintf('Uploading zip... \n');
        uploadFlag = true;
    else
        uploadFlag = false;
    end
    [cloudFolder,zipFileName] =  ...
        sendToCloud(gcp,myScene,'uploadZip',uploadFlag);
    
    % Render (Local)
%     [oi_3d, result] = scene3d.render;
%     ieAddObject(oi_3d);
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
