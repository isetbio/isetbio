function visualizeOTF(theOI, targetWavelength, otfRange, varargin)
p = inputParser;
p.addParameter('axesHandle', [], @ishandle);
p.addParameter('extraData', []);

% Parse input
p.parse(varargin{:});
axesHandle = p.Results.axesHandle;
extraData = p.Results.extraData;

optics = oiGet(theOI, 'optics');
waveOTF = fftshift(opticsGet(optics,'otf data',targetWavelength));
xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');

% 
waveMTF = abs(waveOTF);
[~,idx] = min(abs(ySfCyclesDeg));
mtfSlice = waveMTF(idx,:);


sfTicks = [0.01 0.03 0.1 0.3 1 3 10 30 100];

if (isempty(axesHandle))
    figure(); clf;
    axesHandle = subplot(1,1,1);
    fontSize = 20;
else
    fontSize = 12;
end
axes(axesHandle);

if (~isempty(extraData))
    yyaxis left
end
plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 01.0], 'MarkerSize', 10);
set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
ylabel('modulation');

if (~isempty(extraData))
    yyaxis right
    hold on
    legends = {'PSF'};
    extraDataColors = [1 0 0; 1 0 1; 0 1 0;];
    maxY = 0;
    for dataSetIndex = 1:numel(extraData)
        plot(extraData{dataSetIndex}.sf, extraData{dataSetIndex}.csf/max(extraData{dataSetIndex}.csf), ...
            'rs-', 'MarkerSize', 10, ...
            'MarkerEdgeColor', squeeze(extraDataColors(dataSetIndex,:)),  ...
            'MarkerFaceColor', squeeze(extraDataColors(dataSetIndex,:))*0.5+[0.5 0.5 0.5]);
        legends{numel(legends)+1} = extraData{dataSetIndex}.legend;
    end
    hold off
    set(gca,'YLim', [0 1.05]);
    ylabel(extraData{dataSetIndex}.ylabel);
    legend(legends, 'Location', 'SouthWest');
end

set(gca, 'XLim', [otfRange/100 otfRange], 'XScale', 'log');
set(gca, 'XTick', sfTicks);
axis('square')
grid('on'); box('on');
set(gca, 'FontSize', fontSize);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal'); 

title(sprintf('MTF (%2.0f nm)', targetWavelength));
end

