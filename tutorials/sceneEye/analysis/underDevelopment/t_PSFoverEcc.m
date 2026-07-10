% SkipFile
%% t_PSFoverEcc.m
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
saveDirName = sprintf('eccPSF_%s',currDate);
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


%% Calculate coordinates of each point

ecc = [0 10 20 30]; % desired

% Using Zemax, we convert the above field angles to object height.
% Zemax calculates the angle "with respect to the object space z-axis and
% the paraxial entrance pupil position on the object space z-axis."
% According to the lens prescription, the entrance pupil position is
% located at 3.040 mm from the aperture (toward image space). That
% corresponds to somewhere near the back of the lens.
obj_height = [0 1.763e4 3.640e4 5.774e4]*10^-3; % meters
pointDistance = 100; % In meters(we treat it as essentially inf.)

% We can do a rough sanity check:
field_angle = atand(obj_height./pointDistance);

nRenders = length(obj_height);

%% Find the image height of each point
%{
hLineFig = figure(); hold on;
peakAngle = zeros(1,nRenders);

for ii = 1:nRenders
    
myScene = sceneEye('blackBackdrop');
myScene.name = sprintf('psf_%0.2fdeg',ecc(ii));

myScene.accommodation = 1/pointDistance;
myScene.numCABands = 1;
myScene.pupilDiameter = 4;
myScene.numRays = 128;
myScene.resolution = 4096;

myScene.recipe = recipeSet(myScene.recipe,'cropwindow',[0 1 0.49 0.51]);
pixelFOV = 0.1; 
pointWidth = 2*pointDistance*tand(pixelFOV/2);

% Increase FOV to span all points
full_fov = 2*(max(ecc)+5); 
myScene.fov = full_fov;

myScene.recipe = piAddSphere(myScene.recipe,...
    'spectrum',[400 1 800 1],...
    'radius',pointWidth/2,...
    'location',[-obj_height(ii) pointDistance 0]);

% Render (Local)
[oi, result] = myScene.render;
ieAddObject(oi);
oiWindow;

% Take a horizontal line
photons = oiGet(oi,'photons');
midpt = round(size(photons,1)/2);
wavelengths = oiGet(oi,'wave');
indexGreen = find(wavelengths == 550);
photonsGreen = photons(midpt,1:end,indexGreen);

% Convert to degrees (spherical)
thetaLine = myScene.angularSupport;
figure(hLineFig);
plot(thetaLine,photonsGreen,'k')
grid on; 

% Find peaks for point
peakI = find(photonsGreen == max(photonsGreen));
peakAngle(ii) = thetaLine(peakI);

% Plot peak
plot(peakAngle(ii),max(photonsGreen),'r*');

end
%}
load('peakAngle.mat');
peakAngle

%% Now plot a very zoomed in version for each point
% We can do this more easily now that we have the exact location of each
% point on the image.

for ii = 1:nRenders
    
    % Reload a fresh scene.
    % A black disk located 100 meters away and spanning 120 deg FOV.
    myScene = sceneEye('blackBackdrop');
        
    % Set name
    myScene.name = sprintf('psf_%0.2fdeg',ecc(ii));
    
    % Set fixed parameters
    myScene.accommodation = 1/100;
    myScene.numCABands = 16;
    % myScene.numCABands = 1;
    myScene.pupilDiameter = 4;
    
    % This needs to be high, since the point is small and hard to hit!
    myScene.numRays = 2^14; 
    % myScene.numRays = 1024; % LQ
    
    % Essentially a point. See calculation in s_comparePSF.m
    pixelFOV = 0.0020;
    % pixelFOV = 0.1; % DEBUG
    pointWidth = 2*pointDistance*tand(pixelFOV/2);
    
    % Add equal energy sphere on the optical axis. We will make it the same size as the
    % pixel above.
    myScene.recipe = piAddSphere(myScene.recipe,...
        'spectrum',[400 1 800 1],...
        'radius',pointWidth/2,...
        'location',[-obj_height(ii) pointDistance 0]); % obj_height(ii)
    
    %% Set crop parameters 
    % We do this so we don't have to render the whole retinal image for
    % each point. This could potentially be adapted back into the
    % eccentricity code in sceneEye, which is currently not working.
    
    % We assume we only shift along x eccentricity and y eccentricity is
    % fixed at zero.
    
    % First we calculate the size of the retina that spans all points
    full_fov = 2*(max(ecc)+5); % Add a 5 deg buffer zone
    myScene.fov = full_fov;
    %myScene.fov = 70; %debug
    full_size = myScene.width;

    % Desired FOV for the cropped image
    fov = 0.5;
    % fov = full_fov;
    
    % Convert to image height range
    image_height_x(1) = tand(peakAngle(ii)+fov/2)*myScene.distance2chord;
    image_height_x(2) = tand(peakAngle(ii)-fov/2)*myScene.distance2chord;
    image_height_x = sort(image_height_x); % Low to high
    
    image_height_y(1) = tand(fov/2)*myScene.distance2chord;
    image_height_y(2) = tand(-fov/2)*myScene.distance2chord;
    image_height_y = sort(image_height_y); % Low to high
    
    image_height = [image_height_x image_height_y];
    
    % Convert to image coordinate space
    % (0,0) is upper left hand corner
    image_height  =  image_height + full_size/2;
    
    % Convert from image height range to crop ratio
    crop_window= image_height./full_size;
    myScene.recipe.set('cropwindow',crop_window);
    
    % Calculate resolution
    % Backcalculate the required resolution of the full image, given that
    % we want each crop window to have the folloiwng res:
    psf_res = round(fov/0.0021); % Same as the one in s_comparePSF.m
    % psf_res = 122; % debug
    image_fraction = [crop_window(2)-crop_window(1) crop_window(4)-crop_window(3)];
    full_res = psf_res./image_fraction;
    myScene.resolution = round(full_res(1)); % Assume square
    
    % DEBUG
    % myScene.recipe.set('cropwindow',[0 1 0 1]);
    % myScene.resolution = 256;
    
    %% Render with cloud
    if(ii == nRenders)
        fprintf('Uploading zip... \n');
        uploadFlag = true;
    else
        uploadFlag = false;
    end
    [cloudFolder,zipFileName] =  ...
        sendToCloud(gcp,myScene,'uploadZip',uploadFlag);

    % Render (Local)
%     [oi, result] = myScene.render;
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
