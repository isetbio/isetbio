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

%% Initialize 
clear
% ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

stimulusTestI = 1;
% for stimulusTestI = 1:2
    
switch stimulusTestI
    case 1
        stimulusTest = 'WN';
    case 2 
        stimulusTest = 'NSEM';
end
%% Load stimulus movie using RemoteDataToolbox


switch stimulusTestI
    case 1
        % Binary white noise test movie
%         load(['/Users/james/Documents/MATLAB/'...
%             'akheitman/WN_mapPRJ/Stimuli/'...
%             'BW-8-1-0.48-11111_RNG_16807/testmovie_8pix_Identity_8pix.mat']);

        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');
        testmovie = data.testmovie;
        
%         load('isetbio misc/RDT Uploads/recorded_spikes_fit_WN_test_WN_2013_08_19_6_cell_1203.mat');
%         load('isetbio misc/RDT Uploads/xval_mosaic_WN_ONParasol_2013_08_19_6.mat');
%         load('isetbio misc/RDT Uploads/xval_mosaic_WN_OFFParasol_2013_08_19_6.mat');
        data = rdt.readArtifact('xval_mosaic_WN_ONParasol_2013_08_19_6', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
    case 2
        % NSEM test movie
        % These are small black and white van hatteren images with eye
        % movements superimposed.
        load(['/Users/james/Documents/MATLAB/'...
            'akheitman/NSEM_mapPRJ/Stimuli/'...
            'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);

        % load('isetbio misc/RDT Uploads/recorded_spikes_fit_WN_test_NSEM_2013_08_19_6_cell_1203.mat');
        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_ONParasol_2013_08_19_6.mat');
%         load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_OFFParasol_2013_08_19_6.mat');
end

testmovieshort = testmovie.matrix(:,:,1:1200); 

%% Show test movie
% showFrames = 50;
% ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate outer segment object

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the outer segment responses.  That form of the
% outer segment object is called 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

%% Generate RGC object for simulated GLM prediction of response
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.

% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0; ntrials = 0;

params.experimentID = '2013-08-19-6'; 
params.cellType = 'On Parasol';
params.stimulusTest = stimulusTest; % WN or NSEM, from above

% Create object
innerRetina = irPhys(os1, params);

nTrials = 57;
innerRetina = irSet(innerRetina,'numberTrials',nTrials);

%% Create a new inner retina object and attach the recorded spikes
innerRetinaRecorded = irPhys(os1, params);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

for cellind = 1:length(xval_mosaic)
    for iTrial = 1:nTrials
        recorded_spiketimes{cellind,1,iTrial} = (xval_mosaic{cellind}.rasters.recorded{iTrial});
    end
end

innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',recorded_spiketimes);
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');


%% Compute the inner retina response

% Lienar convolution
innerRetina = irCompute(innerRetina, os1);

% innerRetina = irComputeContinuous(innerRetina, os1);
% % % Spike computation
% for tr = 1:2%ntrials
%     innerRetina = irComputeSpikes(innerRetina, os1);
% end

innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');


%%
% figure; 
% scatter(innerRetinaPSTH{i}(600+(1:minlen-1200)),psth_rec_all{i}(1:minlen-1200),'x');
figure; 
i =2;
minlen = min([length(innerRetinaPSTH{i}) length(innerRetinaRecordedPSTH{i}) ]);
switch stimulusTestI
    case 1
        plot(innerRetinaPSTH{i}(600+(1:minlen-1200)));
    case 2
        plot(innerRetinaPSTH{i}(1200+(1:minlen-1200)));
end
hold on;
% plot(psth_rec_all{i}(1:minlen-1200),'r')
plot(innerRetinaRecordedPSTH{i}(00+(1:minlen-1200)),'r');
%%
for i = 1:2

switch stimulusTestI
    case 1
        rsim = innerRetinaPSTH{i}(600+(1:minlen-1200));
    case 2
        rsim = innerRetinaPSTH{i}(1200+(1:minlen-1200));
end
% rrec = psth_rec_all{i}(1:minlen-1200);
rrec = innerRetinaRecordedPSTH{i}(1:minlen-1200);
Fs(stimulusTestI,i) = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);
% J = ;

end
% end
%%
vcNewGraphWin([],'upperleftbig'); 
hold on;
scatter(Fs(1,:),Fs(2,:),'g','filled');
scatter(Fs(1,:),Fs(2,:),8,'k','filled');

plot(0:.1:1,0:.1:1);
axis([0 1 0 1]);
xlabel('WN Score (AU)'); ylabel('NSEM Score (AU)');
title('Fractional Variance');
set(gca,'fontsize',16);

% 537-(221+127) = 189 total cells used
% 68+118 = 186 cells including bad ones
%%

innerRetinaRecordedOdd = irPhys(os1, params);
innerRetinaRecordedOdd = irSet(innerRetinaRecordedOdd,'numberTrials',nTrials);

iTctr = 0;
for cellind = 1:length(xval_mosaic)
    for iTrial = 1:2:nTrials
        iTctr = iTctr+1;
        recorded_spiketimes{1,cellind,iTctr} = find(xval_mosaic{cellind}.rasters.recorded(iTrial,:)>0);
    end
end

innerRetinaRecordedOdd = irSet(innerRetinaRecordedOdd,'recordedSpikes',recorded_spiketimes);

innerRetinaRecordedOddPSTH = mosaicGet(innerRetinaRecordedOdd.mosaic{1},'responsePsth');


innerRetinaRecordedEven = irPhys(os1, params);
innerRetinaRecordedEven = irSet(innerRetinaRecordedEven,'numberTrials',nTrials);

iTctr = 0;
for cellind = 1:length(xval_mosaic)
    for iTrial = 2:2:nTrials
        iTctr = iTctr+1;
        recorded_spiketimes{1,cellind,iTctr} = find(xval_mosaic{cellind}.rasters.recorded(iTrial,:)>0);
    end
end

innerRetinaRecordedEven = irSet(innerRetinaRecordedEven,'recordedSpikes',recorded_spiketimes);

innerRetinaRecordedEvenPSTH = mosaicGet(innerRetinaRecordedEven.mosaic{1},'responsePsth');

for i = 1:68

J(i) = 1 - sum((innerRetinaRecordedEvenPSTH{i}-innerRetinaRecordedOddPSTH{i}).^2)/sum((innerRetinaRecordedEvenPSTH{i}-mean(innerRetinaRecordedEvenPSTH{i})).^2)


end

% end

%% Plot a few simple properties of the rgcs in the mosaic

% % Spatial RF
% irPlot(innerRetina,'sRFcenter','cell',[1 1]);
% axis([0 13 0 13 -0.1 0.5]); view(0,-90); caxis([-0.2348 0.2348]);
% axis square; axis off; ntrials = 0;
% colormap gray; shading interp
% set(gcf,'position',[ 0.0292    0.5433    0.2951    0.3500])
% % set(gcf,'position',[0.0986    0.6856    0.1944    0.2078]);
% 
% % Temporal impulse response
% irPlot(innerRetina,'tCenter','cell',[1 1]);
% axis square;
% % set(gca,'fontsize',12);
% set(gcf,'position',[    0.3417    0.5422    0.2931    0.3511])
% % set(gcf,'position',[   0.4104    0.6867    0.1556    0.2067])
% 
% % Post spike filter
% irPlot(innerRetina,'postSpikeFilter','cell',[1 1]);
% % set(gca,'fontsize',12);
% axis square;
% set(gcf,'position',[ 0.6590    0.5367    0.2931    0.3544])
% % set(gcf,'position',[   0.6583    0.6856    0.1319    0.2056])

%% Plot the raster and PSTH responses of an RGC 

% Choose cell of interest
i = 2;

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


%% Plot the responses as computed from EJ's lab code
% % Load results locally
% % expdate = '2012-08-09-3'; 
expdate = '2013-08-19-6';
type = 'NSEM';
% glmFitPath = '/Users/james/Documents/matlab/akheitman/NSEM_mapPRJ/';
glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/';
% glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/Test_NSEM/';
matFileNames = dir([glmFitPath expdate '/ON*.mat']);
cell = matFileNames(i).name(1:end-4);
load([glmFitPath expdate '/' cell '.mat']);

% % Load results with RDT
% client = RdtClient('isetbio'); client.crp('/resources/data/rgc');
% [data, artifact] = client.readArtifact('parasol_on_1205', 'type', 'mat');
% fittedGLM = data.fittedGLM;

[psth_sim, psth_rec] = plotrasters(fittedGLM.xvalperformance, fittedGLM);
subplot(211); axis([0 4 0 2*57]); subplot(212); axis([0 4 0 max(psth_sim)]);
 set(gcf,'position',[138          86        1264         500]);