function plotHexMeanImage(conemosaicH, varargin)
% Plots the mean absorptions or current from a coneMosaicHex as a nice
% image. Used in coneMosaicHex.window().
% 
%       plotHexMeanImage(cmosaicH,'type','current')
% 
% This allows coneMosaicWindow to generate the nice images from the
% coneMosaicHex class. This is called from @coneMosaicHex/plot.m which is
% utilized by coneMosaicWindow.m.
% 
% See also t_coneMosaicHex3.m, line 70.
% 
% 10/2016 (c) JRG/BW/NPC isetbio team

%% Parse input
p = inputParser;
addRequired(p,'conemosaicH',@(x) isa(x, 'coneMosaicHex'));
addParameter(p,'type','absorptions',@ischar);
p.parse(conemosaicH, varargin{:});
conemosaicH = p.Results.conemosaicH;
plotType = p.Results.type;

%% Select type of data to plot
switch plotType
    case 'absorptions'
        dataHex = conemosaicH.absorptions;
    otherwise
        dataHex = conemosaicH.current;
end

%% Get the nice coneMosaicHex image for the average mosaic response over time
% Render activation images for the hex mosaic
[activationsHexImage,~,supX,supY] = ...
    conemosaicH.computeActivationDensityMap(mean(dataHex,3));

%% Display results
imagesc(supX,supY,activationsHexImage);
axis 'image'; axis 'xy'
colormap gray;

hold on;
sampledHexMosaicXaxis = conemosaicH.patternSupport(1,:,1) + conemosaicH.center(1);
sampledHexMosaicYaxis = conemosaicH.patternSupport(:,1,2) + conemosaicH.center(2);
dx = conemosaicH.pigment.pdWidth;

axis 'equal'; axis 'xy'
xTicks = [sampledHexMosaicXaxis(1) conemosaicH.center(1) sampledHexMosaicXaxis(end)];
yTicks = [sampledHexMosaicYaxis(1) conemosaicH.center(2) sampledHexMosaicYaxis(end)];
xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], 'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
box on; grid off;
dataRange = prctile(activationsHexImage(:), [5 95]);

set(gca, 'CLim', dataRange);
set(gca, 'XLim', [sampledHexMosaicXaxis(1)-dx sampledHexMosaicXaxis(end)+dx]);
set(gca, 'YLim', [sampledHexMosaicYaxis(1)-dx sampledHexMosaicYaxis(end)+dx]);

drawnow
