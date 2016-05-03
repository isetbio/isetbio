% t_rgcNaturalScenes
% 
% Reproduce Fig. 2 of Heitman, Brackbill, Greschner, Litke, Sher &
% Chichilnisky, 2016 with Isetbio.
% 
% http://biorxiv.org/content/early/2016/03/24/045336
% 
% Load a movie stimulus, generate an os object for the movie, generate an
% inner retina object, load parameters from a physiology experiment in the
% Chichilnisky Lab and compute the response of the RGC mosaic.
% 
% 4/2016
% (c) isetbio team

%% 
clear
% ieInit;


%%
stimulusTestI = 1;
% for stimulusTestI = 1%1:2
    
switch stimulusTestI
    case 1
        stimulusTest = 'WN';
    case 2 
        stimulusTest = 'NSEM';
end
%% Load stimulus movie

% Load stimulus movie using RemoteDataToolbox
% These are small black and white van hatteren images with eye movements
% superimposed.

% rdt = RdtClient('isetbio');
% rdt.crp('resources/data/rgc');
% data = rdt.readArtifact('testmovieshort', 'type', 'mat');
% testmovieshort = data.testmovieshort; 
% % implay(testmovieshort,10);

% Natural scene test movie
% load(['/Users/james/Documents/MATLAB/'...
%     'akheitman/NSEM_mapPRJ/Stimuli/'...
%     'NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat']);

switch stimulusTestI
    case 1
        % White noise test movie
        load(['/Users/james/Documents/MATLAB/'...
            'akheitman/WN_mapPRJ/Stimuli/'...
            'BW-8-1-0.48-11111_RNG_16807/testmovie_8pix_Identity_8pix.mat']);

    case 2
        % NSEM movie
        load(['/Users/james/Documents/MATLAB/'...
            'akheitman/NSEM_mapPRJ/Stimuli/'...
            'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);

%         load(['/Users/james/Documents/MATLAB/'...
%             'akheitman/NSEM_mapPRJ/Stimuli/'...
%             'NSEM_eye-long-v2/X_frame_ONPar_1203.mat']);
%         testmovie.matrix = reshape(X_frame,13,13,601);
end

testmovieshort = testmovie.matrix(:,:,1:1101); 

%% Show test movie
% vcNewGraphWin; 
% for frame1 = 1:200
%     imagesc(testmovieshort(:,:,frame1));
%     colormap gray; 
%     drawnow;
% end
% close;
%% Generate outer segment object

% In this case, the coupled-GLM calculation converts from the frame buffer
% values in the movie to the outer segment responses.  That form of the
% outer segment object is called 'displayRGB'.
os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);
% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

%% Generate RGC object
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.

params.name = 'macaque phys'
params.experimentID = '2013-08-19-6'; % on parasol
% params.outersegment = os1;
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0;

params.cellType = 'On Parasol';
params.stimulusTest = stimulusTest;

innerRetina = irPhys(os1, params);
nTrials = 57;
innerRetina = irSet(innerRetina,'numberTrials',nTrials);


innerRetinaRecorded = irPhys(os1, params);
nTrials = 57;
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

switch stimulusTestI
    case 1
        load('isetbio misc/RDT Uploads/recorded_spikes_fit_WN_test_WN_2013_08_19_6_cell_1203.mat');
    case 2 
        load('isetbio misc/RDT Uploads/recorded_spikes_fit_WN_test_NSEM_2013_08_19_6_cell_1203.mat');
end
innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',recorded_spiketimes);

%% Plot a few simple properties of the rgcs in the mosaic

% Spatial RF
irPlot(innerRetina,'sRFcenter','cell',[1 1]);
axis([0 13 0 13 -0.1 0.5]); view(0,-90); caxis([-0.2348 0.2348]);
axis square; axis off; ntrials = 0;
colormap gray; shading interp
set(gcf,'position',[ 0.0292    0.5433    0.2951    0.3500])
% set(gcf,'position',[0.0986    0.6856    0.1944    0.2078]);

% Temporal impulse response
irPlot(innerRetina,'tCenter','cell',[1 1]);
axis square;
% set(gca,'fontsize',12);
set(gcf,'position',[    0.3417    0.5422    0.2931    0.3511])
% set(gcf,'position',[   0.4104    0.6867    0.1556    0.2067])

% Post spike filter
irPlot(innerRetina,'postSpikeFilter','cell',[1 1]);
% set(gca,'fontsize',12);
axis square;
set(gcf,'position',[ 0.6590    0.5367    0.2931    0.3544])
% set(gcf,'position',[   0.6583    0.6856    0.1319    0.2056])
%% Compute the inner retina response

% Lienar convolution
innerRetina = irCompute(innerRetina, os1);

% % Spike computation
for tr = 1:ntrials
    innerRetina = irComputeSpikes(innerRetina, os1);
end

%% Plot the raster and PSTH responses of an RGC 

% Choose cell of itnerest
i = 1;


irPlot(innerRetinaRecorded ,'raster','cell',[i 1],'color','k');
axis([0 8 0 nTrials]);
% set(gcf,'position',[  0.0965    0.6144    0.8757    0.2622]);
set(gcf,'position',[  0.0896    0.3800     0.8764    0.2633]);
title(sprintf('Recorded Spikes %s Type %s Fit %s Test %s Cell %s', ...
    innerRetina.mosaic{1}.experimentID, innerRetina.mosaic{1}.cellType, innerRetina.mosaic{1}.stimulusFit,...
    innerRetina.mosaic{1}.stimulusTest, strrep(innerRetina.mosaic{1}.cellID{i},'_','\_')));
if strcmpi(stimulusTest,'WN'); axis([0.5 8.5 0 57]); end;

irPlot(innerRetina,'raster','cell',[i 1]);
axis([1 9 0 nTrials]);
% set(gcf,'position',[  0.0965    0.6144    0.8757    0.2622]);
set(gcf,'position',[ 0.0896    0.0733    0.8771    0.2644]);
title(sprintf('Simulated Spikes %s Type %s Fit %s Test %s Cell %s', ...
    innerRetina.mosaic{1}.experimentID, innerRetina.mosaic{1}.cellType, innerRetina.mosaic{1}.stimulusFit,...
    innerRetina.mosaic{1}.stimulusTest, strrep(innerRetina.mosaic{1}.cellID{i},'_','\_')));

%% Plot the PSTH responses for RGCs

% Load precomputed results from EJ's lab code
% load('isetbio misc/scratch/xvalall_59_trials2.mat');
% load('isetbio misc/scratch/psth_rec_all.mat');

% Load simulated PSTHs from EJ's code
% rdt = RdtClient('isetbio'); rdt.crp('resources/data/rgc');
% data = rdt.readArtifact('xvalall_59_trials2', 'type', 'mat');
% xvalall = data.xvalall;
% 
% % Load recorded PSTHs from EJ's code
% rdt = RdtClient('isetbio'); rdt.crp('resources/data/rgc');
% data = rdt.readArtifact('psth_rec_all', 'type', 'mat');
% psth_rec_all = data.psth_rec_all;

switch stimulusTestI
    case 1
        load('isetbio misc/RDT uploads/WN_psth_lab_rec_2013_08_19_6_cell_1203_long.mat');
    case 2        
        load('isetbio misc/RDT uploads/WN_Test_NSEM_psth_lab_rec_2013_08_19_6_cell_1203_long.mat');
end

% load('isetbio misc/RDT uploads/psth_lab_rec_2013_08_19_6_cell_1203.mat');

psth_rec_all{i} = xvalall{2}.psth_rec;
% subplot(212);
vcNewGraphWin;
rgc2psth = mosaicGet(innerRetina.mosaic{1},'responsePsth');

% % rgc2linear = mosaicGet(innerRetina.mosaic{1},'responseLinear');
% rgc2linear = innerRetina.mosaic{1}.responseLinear;
% figure; plot(rgc2linear{1,1,1})

% minlen = min([length(xvalall{2}.psth)]);
minlen = min([ length(psth_rec_all{i}) length(rgc2psth{i})]);
    
hold on;

% Plot isetbio calculation
plot([1:minlen-1200]./1208,rgc2psth{i}(1200+(1:minlen-1200)),'r ','linewidth',3);
% irPlot(innerRetina,'psth','cell',[i 1]); hold on;

% Plot output of Chichilnisky Lab code
plot((1200+[1:minlen-1200])./1208,1*xvalall{2}.psth(1200+(1:minlen-1200)),':b','linewidth',2);
% % plot([1:minlen-1200]./1208,(20/59)*xvalall{i}.psth(00+(1:minlen-1200)),':k','linewidth',2);

% Plot recorded response to natural scene
plot([1:minlen-1200]./1208,(1)*psth_rec_all{i}(1:minlen-1200),'k','linewidth',2);

% Set display properties
% set(gcf,'position',[0.0931    0.2856    0.8806    0.2533]);
set(gcf,'position',[ 0.0931    0.0300    0.8764    0.2633]);
% axis([0 (length(rgc2psth{i}))./1208 0  max(rgc2psth{i})]);%max(rgc2psth{i}(1:minlen))])
axis([0 8 0  max(rgc2psth{i})]);%max(rgc2psth{i}(1:minlen))])
[maxv, maxi] = max(rgc2psth{i}(1:minlen)-xvalall{2}.psth(1:minlen));
title(sprintf('Validation with Chichilnisky Lab Data\nCell %d, maxv = %.1f, maxi = %d',i,maxv,maxi));
xlabel('Time (sec)'); ylabel('PSTH (spikes/sec)');
legend('ISETBIO','Lab Code','Recorded');
% legend('ISETBIO','Recorded');
set(gca,'fontsize',14)



%%

% end


%%
% figure; 
% scatter(rgc2psth{i}(600+(1:minlen-1200)),psth_rec_all{i}(1:minlen-1200),'x');
figure; 
plot(rgc2psth{i}(600+(1:minlen-1200)));
hold on;
plot(psth_rec_all{i}(1:minlen-1200),'r')

rsim = rgc2psth{i}(600+(1:minlen-1200));
rrec = psth_rec_all{i}(1:minlen-1200);
Fs = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);
% J = ;
%% Plot the responses as computed from EJ's lab code
% % Load results locally
% % expdate = '2012-08-09-3'; 
expdate = '2013-08-19-6';
type = 'NSEM';
% glmFitPath = '/Users/james/Documents/matlab/akheitman/NSEM_mapPRJ/';
glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/';
% glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/Test_NSEM/';
matFileNames = dir([glmFitPath expdate '/ON*.mat']);
cell = matFileNames(2).name(1:end-4);
load([glmFitPath expdate '/' cell '.mat']);

% % Load results with RDT
% client = RdtClient('isetbio'); client.crp('/resources/data/rgc');
% [data, artifact] = client.readArtifact('parasol_on_1205', 'type', 'mat');
% fittedGLM = data.fittedGLM;

[psth_sim, psth_rec] = plotrasters(fittedGLM.xvalperformance, fittedGLM);
subplot(211); axis([0 4 0 2*57]); subplot(212); axis([0 4 0 max(psth_sim)]);
 set(gcf,'position',[138          86        1264         500]);