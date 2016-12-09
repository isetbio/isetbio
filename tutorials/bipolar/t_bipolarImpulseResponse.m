%% t_bipolarImpulseResponse
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses.
% 
% The stimulus is an impulse. The goal is to examine the impulse repsonse
% of the bipolar object with the differentiator 
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
ieInit

%% Create impulse ois





% Set photon rates. This is a kluge that appeared
% just for this test, and that should probably go
% away again. This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Create outersegment object and get the adapted response.
noiseFlag = 0;
if strcmp(osModelType,'linear')
    os = osLinear(sensor); 
    paramsOS.convolutionType = 0; 
else
    os = osBioPhys(sensor); 
    paramsOS.bgVolts = 10*mean(vectorize(sensorGet(sensor,'volts')));
end

os = osSet(os, 'noiseFlag', noiseFlag);
os = osCompute(os, sensor, paramsOS);

% Set size of retinal patch
patchSize = sensorGet(sensor,'width','m');
os = osSet(os, 'patch size', patchSize);

% % Set time step of simulation equal to absorptions
% timeStep = sensorGet(absorptions,'time interval','sec');
% os = osSet(os, 'time step', timeStep);
        
% vcNewGraphWin; plot(squeeze(stimulus))
% xlabel('Time (msec)','fontsize',14); ylabel('Stimulus Intensity','fontsize',14)
% osPlot(os,sensor);

% Plot all cone responses
% ccrs = reshape(os.coneCurrentSignal,64*64,400);
% figure; plot(ccrs')

%% Find bipolar responses


bpParams.filterType = 4;
bpParams.cellLocation = 58;
bp = bipolar(os,bpParams);
bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);
% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, os);

% bipolarPlot(bp,'response');

bpResponse = bipolarGet(bp,'response');

tInd = 100;
bpRFrame = bpResponse(:,:,tInd);
[maxbpR,maxbpRind] = max(bpRFrame(:));
[rmax,cmax] = ind2sub(size(bpRFrame),maxbpRind);

tBin = osGet(os, 'timeStep');
vcNewGraphWin([],'upperleftbig'); 
plot(tBin/.001:tBin/.001:(tBin/.001)*size(bpResponse,3),squeeze(bpResponse(rmax,cmax,:) - bpResponse(rmax,cmax,end))./max(abs(squeeze(bpResponse(rmax,cmax,:) - bpResponse(rmax,cmax,end)))));
% figure; plot(tBin/.001:tBin/.001:(tBin/.001)*length(bpResponse),squeeze(bpResponse - bpResponse(end)));

load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/OnParasolExcFilters.mat')
% figure; 
hold on;
% plot(tme,max(abs(1))*mean(LinFilt)./max(abs(mean(LinFilt))));

meanLinFilt = mean(LinFilt);
% Downsample to appropriate time course
dsrate = length(tme)/max(tme);
rgcFilt = meanLinFilt(1:2*dsrate:end);
% plot(tme(1:dsrate:length(meanLinFilt)),meanLinFilt(1:dsrate:end));
bpComp = squeeze(bpResponse(rmax,cmax,:) - bpResponse(rmax,cmax,end))./max(abs(squeeze(bpResponse(rmax,cmax,:) - bpResponse(rmax,cmax,end))));
plot(tme,meanLinFilt);
xlabel('Time (msec)','fontsize',16);
ylabel('Normalized Amplitude','fontsize',16);

% figure; plot((bpComp(1:1*dsrate:1*1880))); hold on; plot(18+(1:189),rgcFilt(1:189))
% corrcoef((bpComp(1:1*dsrate:1*1880)),[zeros(1,18) rgcFilt(1:189-19)]')
dsrate2 = .001/timeStep;
% figure; plot(bpComp(1:2*dsrate2:2*188*dsrate2)); hold on; plot(rgcFilt(1:188));
cM = corrcoef((bpComp(1:2*dsrate2:2*188*dsrate2)),rgcFilt(1:188));

title(sprintf('On Parsol Synaptic Input Impulse Response\nR^2 = %0.2f',cM(2,1)),'fontsize',16);
legend('Bipolar IR','RGC Synaptic IR');
set(gca,'fontsize',16);
grid on;


%%


%% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaSU = irCreate(bp, params);

% innerRetinaSU.mosaicCreate('model','glm','type','off midget');
innerRetinaSU.mosaicCreate('model','glm','type','on midget');

innerRetinaSU = irCompute(innerRetinaSU, bp);
irPlot(innerRetinaSU,'linear');
%%
% irPlot(innerRetinaSU,'linear')
% %% Find RGC responses
% % Build and IR object that takes as input the bipolar mosaic.
% 
% % Initialize.
% clear params 
% clear innerRetinaBpSu
% params.name      = 'Bipolar with nonlinear subunits'; % This instance
% params.eyeSide   = 'left';   % Which eye
% params.eyeRadius = 4;        % Radius in mm
% params.eyeAngle  = 90;       % Polar angle in degrees
% % bp = bipolarSet(bp,'patchSize',2e-4);
% innerRetinaBpSu = irCreate(bp, params);
% 
% % Create a subunit model for the on midget ganglion cell parameters
% innerRetinaBpSu.mosaicCreate('model','Subunit','type','off parasol');
% innerRetinaBpSu.mosaic{1}.mosaicSet('numberTrials',10);
% % % Uncomment to get rid of spatial nonlinearity
% % newRectifyFunction = @(x) x;
% % innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);
% 
% innerRetinaBpSu.mosaic{1}.mosaicSet('tonicDrive',0.01);
% 
% % irPlot(innerRetinaBpSu,'mosaic');
% 
% % Compute RGC mosaic responses
% innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
% irPlot(innerRetinaBpSu,'linear')

%% Measure the temporal impulse response for the original GLM model
% This response can be compared to the above response 
% Image parameters
params.image_size = 63;
params.meanLuminance = 1;
params.nsteps = 100;
params.fov = 0.8;

% Set impulse stimulus by only turning on a signal pixel at t = 1
sceneRGB = zeros(params.image_size, params.image_size, params.nsteps, 3);
sceneRGB(16,16,23,:) = 100*[10 10 10];

% Create outer segment
osD = osCreate('displayRGB');

% Set os parameters
coneSpacing = 100e-6;
osD = osSet(osD, 'patchSize', coneSpacing);

coneSampling = 0.002;
osD = osSet(osD, 'timeStep', coneSampling);

% Set os stimulus
osD = osSet(osD, 'rgbData', sceneRGB);

% Build inner retina
clear params innerRetinaRGB
params.name      = 'Macaque inner retina 1 impulse'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaRGB = irCreate(osD, params);

innerRetinaRGB.mosaicCreate('model','glm','type','off parasol');

% Compute response
innerRetinaRGB = irCompute(innerRetinaRGB, osD);

% Plot response
irPlot(innerRetinaRGB,'linear');%,'cell',[1 2])
axis([0 0.45 2.245 2.28])

%%
os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', sceneRGB);
