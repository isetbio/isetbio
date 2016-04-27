% v_rgcExternal
% 
% A validation for the coupled GLM model of retinal ganglion cells.
% This file compares finds the response of an ON Parsol RGC mosaic using
% both the isetbio implementation of the GLM and the Chichilnisky Lab's
% code. The ON Parasol RGC parameters were fit to physiological data from
% the Chichilnisky lab labeled as 2012-08-09-3.
% 
% The Remote Data Toolbox is used to download the necessary data from an
% Archiva server. Access is password protected; email the isetbio contact
% for access. Make sure you have the Remote Data Toolbox repository 
% (available on the isetbio github at isetbio/RemoteDataToolbox) added to
% your Matlab path. Four files are loaded using the RDT:
% 
%       1. testmovishort: The movie stimulus. In this case, a short 
%       sequence of natural scenes with simulated eye movements.
%       2. parasol_on_1205: The fitted parameters for the ON Parasol mosaic. 
%       This consists  of 39 separate "fittedGLM" structures.
%       ******** Right now, this is only 1 cell; will change soon.*******
%       3. pairspikesall: The measured spiking activity of the mosaic. This 
%       input is required for the coupled GLM.
%       4. xvalall: The output of the same simulation from code used by the
%       Chichilnisky Lab.
% 
% The script uses the osIdentity object for the cone outersegment
% representation, which simply copies the RGB data from the isetbio scene
% for each frame of the movie. This is due to the fact that the GLM
% is input-referred as opposed to cone current-referred. A new version of
% the isetbio GLM implementation will be cone current-referred and will be
% added in the near future.
% 
% (c) isetbio team
% JRG 12/2015
%% 
clear
ieInit;

%% Load stimulus movie

% Load stimulus movie using RemoteDataToolbox
% These are small black and white van hatteren images with eye movements
% superimposed.
% rdt = RdtClient('isetbio');
% rdt.crp('resources/data/rgc');
% data = rdt.readArtifact('testmovieshort', 'type', 'mat');
% testmovieshort = data.testmovieshort;
% % implay(testmovieshort,10);

% load('/Users/james/Documents/MATLAB/akheitman/NSEM_mapPRJ/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat');
load('isetbio misc/scratch/testmovieshort.mat');
% testmovieshort = testmovie.matrix(:,:,1:601);


% figure;
% for frame1 = 1:200
%     imagesc(testmovieshort(:,:,frame1));
%     colormap gray; 
%     drawnow;
% end
% close;
%% Generate outer segment object

% In this case, the coupled-GLM calculation converts from the frame buffer
% values in the movie to the outer segment responses.  That form of the
% outer segment object is called 'identity'.
os2 = osCreate('displayrgb');

% Attach the move to the object
os2 = osSet(os2, 'rgbData', double(testmovieshort));

%% Generate RGC object
params.name = 'macaque phys'
params.outersegment = os2;
params.eyeSide = 'left'; 
params.eyeRadius = 9; 
params.eyeAngle = 90;
% rgc2 = rgcCreate('rgcPhys', params);
rgc2 = irPhys(os2, params);
% remove number of trials!

nTrials = 59;
rgc2 = irSet(rgc2,'numberTrials',nTrials);

%%
rgc2 = irCompute(rgc2, os2);

% for tr = 1:nTrials
%     rgc2 = irComputeSpikes(rgc2, os2);
% end

% rgc2Linear = mosaicGet(rgc2.mosaic{1},'responseLinear');
rgc2psth = mosaicGet(rgc2.mosaic{1},'responsePsth');

%% Load validation data
% % % Load RDT version of output from the Chichilnisky Lab's code
% rdt = RdtClient('isetbio');
% % client.credentialsDialog();
% rdt.crp('resources/data/rgc');
% [data, artifact] = rdt.readArtifact('xvalall', 'type', 'mat');
% xvalall = data.xvalall;

load('isetbio misc/scratch/xvalall_59_trials2.mat');

% load('isetbio misc/scratch/xvalall_20_trials.mat');

% % Load local copy
% load('xvalall2.mat');
load('isetbio misc/scratch/psth_rec_all.mat');
% load('isetbio misc/scratch/psth_rec_all_59_trials.mat');
% load('isetbio misc/scratch/psth_rec_all_20_trials.mat');
%% Compare isetbio output and Chichilnisky Lab output
vcNewGraphWin([],'upperleftbig');
for i = 1%:36
    % Measure difference between outputs
    minlen = min([length(rgc2psth{i}) length(xvalall{i}.psth)]);
    diffpsth(i) = sum(abs(rgc2psth{i}(1:minlen) - xvalall{i}.psth(1:minlen)))./sum(.5*(rgc2psth{i}(1:minlen) + xvalall{i}.psth(1:minlen)));
    
    % % Plot difference
    % subplot(6,7,i); hold on;
    % plot(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen),'b','linewidth',1);

    hold on;
    % Plot recorded response to natural scene
    plot([1:minlen-1200]./1208,(20/59)*psth_rec_all{i}(1:end-1200),'b','linewidth',2);
    
    % Plot output of isetbio code
    % plot([1:minlen]./1208,rgc2psth{i}((1:minlen)),'r ','linewidth',3);
    plot([1:minlen-1200]./1208,rgc2psth{i}(1200+(1:minlen-1200)),'r ','linewidth',3);
    
    % Plot output of Chichilnisky Lab code
    plot([1:minlen-1200]./1208,xvalall{i}.psth(1200+(1:minlen-1200)),':k','linewidth',2);
%     plot([1:minlen-1200]./1208,(20/59)*xvalall{i}.psth(00+(1:minlen-1200)),':k','linewidth',2);

    
    axis([0 (6285-1200)./1208 0  1500]);%max(rgc2psth{i}(1:minlen))])
    [maxv, maxi] = max(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen)); 
    title(sprintf('Validation with Chichilnisky Lab Data\nCell %d, maxv = %.1f, maxi = %d',i,maxv,maxi));
    xlabel('Time (sec)'); ylabel('PSTH (spikes/sec)');
%     legend('Recorded','ISETBIO','Lab Code');
%     legend('ISETBIO','Lab Code');
    set(gca,'fontsize',14)
end

% irPlot(rgc2, 'psth','cell',[1 1]);

% irPlot(rgc2, 'raster','cell',[1 1]);

% figure; plot(horzcat(rgc2.mosaic{1}.spikeResponse{:,1,1,3}))

% vcNewGraphWin; plot(diffpsth,'x')
%%
% figure; % set(gcf,'position',[0    0.0711    1.0000    0.8222]);
vcNewGraphWin([],'upperleftbig');
for i = 1:length(rgc2.mosaic{1}.cellLocation)
    % Measure difference between outputs
    minlen = min([length(rgc2psth{i}) length(xvalall{i}.psth)]);
    diffpsth(i) = sum(abs(rgc2psth{i}(1:minlen) - xvalall{i}.psth(1:minlen)))./sum(.5*(rgc2psth{i}(1:minlen) + xvalall{i}.psth(1:minlen)));
    
%     % % Plot difference
    subplot(6,7,i); hold on;
%     % plot(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen),'b','linewidth',1);
% 
%     % Plot output of isetbio code
%     plot([1:minlen]./1208,rgc2psth{i}(1:minlen),'r ','linewidth',3);
%     hold on;
%     % Plot output of Chichilnisky Lab code
%     plot([1:minlen]./1208,xvalall{i}.psth(1:minlen),':k','linewidth',2);
    hold on;
    % Plot recorded response to natural scene
    plot([1:minlen-1200]./1208,(20/57)*psth_rec_all{i}(1:end-1200),'b','linewidth',2);
    
    % Plot output of isetbio code
    % plot([1:minlen]./1208,rgc2psth{i}((1:minlen)),'r ','linewidth',3);
    plot([1:minlen-1200]./1208,rgc2psth{i}(1200+(1:minlen-1200)),'r ','linewidth',3);
    
    % Plot output of Chichilnisky Lab code
    plot([1:minlen-1200]./1208,xvalall{i}.psth(1200+(1:minlen-1200)),':k','linewidth',2);

    

%     axis([0 6285./1208 0  100]);%max(rgc2psth{i}(1:minlen))])
    [maxv, maxi] = max(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen)); 
    title(sprintf('maxv = %.1f, maxi = %d',maxv,maxi));
%     xlabel('Time (sec)'); ylabel('PSTH (spikes/sec)');
%     legend('ISETBIO','Lab');
%     set(gca,'fontsize',14)
end
%% Compare linear outputs
figure;
for i = 1%:39
%     subplot(6,7,i);
%     plot(rgc2psth{i});
%     hold on;
%     plot(xvalall{i}.psth);
    
    plot(rgc2linear{i}(1:minlen/10)-xvalall{i}.lcif_const(1:10:minlen),'b','linewidth',2);
    hold on;
%     plot(xvalall{i}.cif0);
    
end
