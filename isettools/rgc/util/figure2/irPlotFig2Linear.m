function irPlotFig2Linear(innerRetina, cellInd)

%% Plot a few simple properties of the rgcs in the mosaic

if strcmpi(ieParamFormat(innerRetina.mosaic{1}.cellType),'offparasol')
    offMult = -1; offSwitch = 1;
else
    offSwitch = 0;
end

% Spatial RF
irPlot(innerRetina,'sRFcenter','cell',[cellInd 1]);
if offSwitch
    axis([0 13 0 13 -0.5 0.1]); view(0,-90); caxis([-0.2348 0.2348]);
else 
    axis([0 13 0 13 -0.1 0.5]); view(0,-90); caxis([-0.2348 0.2348]);
end
axis square; axis off; ntrials = 0;
colormap gray; shading interp
set(gcf,'position',[ 0.0292    0.5433    0.2951    0.3500])
% set(gcf,'position',[0.0986    0.6856    0.1944    0.2078]);

% Temporal impulse response
irPlot(innerRetina,'tCenter','cell',[cellInd 1]);
% if offSwitch; view(0,-90); end;
if offSwitch
    axis([0    0.2500   -1.4000    0.4000]);
else
    axis([0    0.2500   -0.4000    1.4000]);
end
axis square;
% set(gca,'fontsize',12);
set(gcf,'position',[    0.3417    0.5422    0.2931    0.3511])
% set(gcf,'position',[   0.4104    0.6867    0.1556    0.2067])

% Post spike filter
irPlot(innerRetina,'postSpikeFilter','cell',[cellInd 1]);
% set(gca,'fontsize',12);
axis square;
set(gcf,'position',[ 0.6590    0.5367    0.2931    0.3544])
% set(gcf,'position',[   0.6583    0.6856    0.1319    0.2056])

drawnow;