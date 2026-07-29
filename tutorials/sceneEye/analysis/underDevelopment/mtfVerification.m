%% mtfVerification.m
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;
if ~mcDockerExists, mcDockerConfig; end % check whether we can use docker
if ~mcGcloudExists, mcGcloudConfig; end % check whether we can use google cloud sdk;

%% Initialize save folder
% Since rendering these images often takes a while, we will save out the
% optical images into a folder for later processing.
currDate = datestr(now,'mm-dd-yy_HH_MM');
saveDirName = sprintf('mtfVer_%s',currDate);
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

clusterName = 'trisha';
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

%% Render a fast image of the slanted bar first

% The slanted bar scene consists of a square plane (1x1 m) that is
% split in half diagonally. The bottom left half is white while the top
% right half is black. By default the plane is placed at [0 0 1] meters,
% but we can change that by given sceneEye an optional 'planeDistance'
% input. 
myScene = sceneEye('slantedBar','planeDistance',50); 
myScene.name = 'mtfVerification';
myScene.numRays = 2048;
myScene.resolution = 512; 

myScene.numBounces = 1;
myScene.numCABands = 16;

myScene.accommodation = 0;
myScene.pupilDiameter = 4;
myScene.fov = 0.5;

[cloudFolder,zipFileName] =  ...
    sendToCloud(gcp,myScene,'uploadZip',true);

% vcAddObject(oi);
% oiWindow;

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

%% Calculate MTF

close all

wavelengths = [400 450 550 650 680 700];

for ii = 1:length(wavelengths)

[freq,mtf] = calculateMTFeye(oi,wavelengths(ii));

plot(freq,mtf); hold on;

end

legendCell = cellstr(num2str(wavelengths', '%d'));
legend(legendCell);

xlabel('Spatial Frequency (cycles/deg)');
ylabel('Contrast Reduction (SFR)');
grid on;
axis([0 60 0 1])
title('MTF (632 nm, 4 mm, 0 dpt)')

set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
