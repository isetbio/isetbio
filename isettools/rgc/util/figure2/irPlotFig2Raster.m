function irPlotFig2Raster(innerRetina, innerRetinaRecorded, cellind, stimulusTestI);
% Plot the rasters for Fig 2.
% 
% Inputs: simulated and recorded IR objects
% 
% 5/2016 JRG (c) isetbio team

% Choose cell of interest
i = cellind;

% Plot recorded spikes
irPlot(innerRetinaRecorded ,'raster','cell',[i 1],'color','k');
nTrials = innerRetina.mosaic{1}.numberTrials;
axis([0 8 0 nTrials]);
% set(gcf,'position',[  0.0965    0.6144    0.8757    0.2622]);
set(gcf,'position',[  0.0896    0.3800     0.8764    0.2633]);
title(sprintf('Recorded Spikes %s Type %s Fit %s Test %s Cell %s', ...
    innerRetina.mosaic{1}.experimentID, innerRetina.mosaic{1}.cellType, innerRetina.mosaic{1}.stimulusFit,...
    innerRetina.mosaic{1}.stimulusTest, strrep(innerRetina.mosaic{1}.cellID{i},'_','\_')));
if stimulusTestI==1; axis([0.5 8.5 0 57]); end;
axis off;

% Plot simulated spikes
set(gca,'fontsize',16');
irPlot(innerRetina,'raster','cell',[i 1]);
axis([1 9 0 nTrials]);
% set(gcf,'position',[  0.0965    0.6144    0.8757    0.2622]);
set(gcf,'position',[ 0.0896    0.0733    0.8771    0.2644]);
title(sprintf('Simulated Spikes %s Type %s Fit %s Test %s Cell %s', ...
    innerRetina.mosaic{1}.experimentID, innerRetina.mosaic{1}.cellType, innerRetina.mosaic{1}.stimulusFit,...
    innerRetina.mosaic{1}.stimulusTest, strrep(innerRetina.mosaic{1}.cellID{i},'_','\_')));
set(gca,'fontsize',16');
axis off

drawnow;