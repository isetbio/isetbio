function filterOut = bipolarFilter(obj, cmosaic, varargin)
% 
% Builds the bipolar temporal filter by deconvolving the outer segment
% impulse response from the RGC impulse response.
% 
% Called internally from bipolarCompute.m
% 
%   filterOut = bipolarFilter(bp, cmosaic);
% 
% (c) isetbio team JRG/BW 12/2016

%% parse input parameters
p = inputParser;
p.addRequired('obj', @(x) (isa(x, 'bipolar')));
p.addRequired('cmosaic', @(x) (isa(x, 'coneMosaic')));  
p.parse(obj, cmosaic, varargin{:});

%% Get cone linear filters

lmsFilters = cmosaic.os.lmsConeFilter;

%% Load RGC object to get time course

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.
testmovieshort = rand(80,40,1);

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep',cmosaic.os.timeStep);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

% Generate RGC object for simulated GLM prediction of response
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.

% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = sqrt(sum(0.^2)); 
params.eyeAngle = 0; ntrials = 0;

% % Determined at beginning to allow looping
% params.experimentID = '2013-08-19-6'; % Experimental dataset
% params.stimulusTest = 'WN'; % WN or NSEM
% params.cellType = 'On Parasol';         % ON or OFF Parasol
% innerRetina = irPhys(os1, params);
% dsrate = round((1/150)/cmosaic.os.timeStep);
% warning('time step must be integer multiple');
% osFilt = lmsFilters(1:dsrate:end,1);

cellType = 'onParasol';
% cellType = 'offParasol';
innerRetina = ir(obj, params);
innerRetina.mosaicCreate('type',cellType,'model','GLM');
osFilt = lmsFilters(:,1);

%% Deconvolve in Fourier domain to compute bipolar IR

% bipolarFiltMat = zeros(length(innerRetina.mosaic{1}.sRFcenter),89);
% for cellNum = 1%:length(innerRetina.mosaic{1}.sRFcenter)

cellNum = 1
rgcFilt = innerRetina.mosaic{1}.tCenter{cellNum};
rgcFiltPad = [rgcFilt; zeros(length(osFilt(1:end))-1,1)];
osFiltPad = [osFilt; zeros(length(rgcFilt)-1,1)];
% figure; plot(osFiltPad); hold on; plot(rgcFiltPad)

% Convolve with a Gaussian window for low-pass filter
% Bigger variance keeps high freq
gaussVar = 200;
gw1 = gausswin(length(rgcFiltPad),gaussVar);
% gw2 = gausswin(length(rgcFiltPad),gaussVar);

% % No lowpass filter
% bipolarFilt = ifft(fft(rgcFiltPad)./(fft(osFiltPad)));

% % Lowpass only RGC Synaptic IR
bipolarFilt = ifft(fft(fftshift(gw1./norm(gw1)^2)).*fft(rgcFiltPad)./fft(osFiltPad));
% bipolarFilt = ifft(fft(rgcFiltPad)./fft(osFiltPad));

% bipolarFiltPad = [bipolarFilt; zeros(length(osFilt(1:end))-1,1)];
rgcOut = ifft(fft(bipolarFilt).*fft(osFiltPad));

%     bipolarFiltMat(cellNum,:) = bipolarFilt;

% end
bipolarFiltMat = bipolarFilt;
% figure; plot(cmosaic.os.timeStep*(1:size(bipolarFiltMat,1)),(bipolarFiltMat)')

%% Compare rgc IR with IR from product of os IR and bipolar IR
dt = cmosaic.os.timeStep;
figure;  plot(dt*(1:size(rgcOut,1)),ieScale(rgcOut),'linewidth',4); 
hold on; plot(dt*(1:size(rgcFiltPad,1)),ieScale(rgcFiltPad),':r','linewidth',3);
legend('os IR * bipolar IR', 'RGC IR');
xlabel('Time (sec)'); ylabel('Normalized Response');
set(gca,'fontsize',14); grid on;

%% Plot cone, bipolar and RGC IRs together

figure; 
% osSigNorm = (squeeze(lmsFilters(:,1)))./max((squeeze(lmsFilters(:,1)) - squeeze(lmsFilters(:,1))));
osSigNorm = ieScale(lmsFilters(:,1));
plot(dt*(1:size(osSigNorm,1)),osSigNorm,'r','linewidth',3)

hold on; plot(dt*(1:length(rgcFiltPad)),rgcFiltPad./max(rgcFiltPad),'g','linewidth',3);
plot(cmosaic.os.timeStep*(0:size(bipolarFiltMat,1)),[0 (bipolarFiltMat)']./max(abs((bipolarFiltMat)')),'b','linewidth',3)

legend('OS IR','RGC IR','Bipolar IR');
title('Bipolar impulse response from deconvolution');
axis([0 0.4 -0.4 1]); grid on;
xlabel('Time (sec)'); ylabel('Normalized Response');
set(gca,'fontsize',14);

filterOut = bipolarFiltMat;