% get data from a figure

function [x,y] = rgcGetDataFromFigure(figureName, closeFigure)
% This returns the contrast values used in the figure (the x values) and
% the probability of correct guess (the y value).

%figureName = '~/RA/results/spikeRateGainTuning/contrastSensitivity/contrastSensitivity_cpd10_SRG80_fullContrastScale.fig';
if notDefined('closeFigure')
    closeFigure = 1;
end

fighandle=openfig(figureName);
ax=findall(fighandle,'Type','line');
x=get(ax,'Xdata');
y=get(ax,'YData');

if closeFigure
    close(fighandle)
end

end