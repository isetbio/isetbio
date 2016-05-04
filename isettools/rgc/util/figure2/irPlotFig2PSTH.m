function irPlotFig2PSTH(innerRetina, psthSim, psthRecorded, cellind, stimulusTestI)
% Plot the PSTHs for Fig 2.
% 
% Inputs: simulated and recorded PSTHs
% 
% 5/2016 JRG (c) isetbio team

vcNewGraphWin([],'upperleftbig'); 

% Choose cell to plot
i = cellind;

minlen = min([length(psthSim{i}) length(psthRecorded{i}) ]);

switch stimulusTestI
    case 1
        plot((00+[1:minlen-1200])./1208, psthSim{i}(600+(1:minlen-1200)),'r','linewidth',3);
    case 2
        plot((00+[1:minlen-1200])./1208, psthSim{i}(1200+(1:minlen-1200)),'r','linewidth',3);
end

hold on;
% plot(psth_rec_all{i}(1:minlen-1200),'r')
plot([1:minlen-1200]./1208,psthRecorded{i}((1:minlen-1200)),':k','linewidth',2);


% Set display properties
% set(gcf,'position',[0.0931    0.2856    0.8806    0.2533]);
set(gcf,'position',[ 0.0931    0.0300    0.8764    0.2633]);
% axis([0 (length(rgc2psth{i}))./1208 0  max(rgc2psth{i})]);%max(rgc2psth{i}(1:minlen))])
axis([0 8 0  max(psthSim{i})]);%max(rgc2psth{i}(1:minlen))])
[maxv, maxi] = max(psthSim{i}(1:minlen));
% title(sprintf('Validation with Chichilnisky Lab Data\nCell %d, maxv = %.1f, maxi = %d',i,maxv,maxi));
title(sprintf('PSTH %s %s, %s Fit, %s Test, Cell %s', ...
    innerRetina.mosaic{1}.experimentID, innerRetina.mosaic{1}.cellType, innerRetina.mosaic{1}.stimulusFit,...
    innerRetina.mosaic{1}.stimulusTest, strrep(innerRetina.mosaic{1}.cellID{i},'_','\_')));
xlabel('Time (sec)'); ylabel('PSTH (spikes/sec)');
% legend('ISETBIO','Lab Code','Recorded');
legend('ISETBIO','Recorded');
set(gca,'fontsize',14)