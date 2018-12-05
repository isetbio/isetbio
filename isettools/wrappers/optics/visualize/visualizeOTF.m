function visualizeOTF(theOI, targetWavelength, otfRange, varargin)
p = inputParser;
p.addParameter('axesHandle', [], @ishandle);
% Parse input
p.parse(varargin{:});
axesHandle = p.Results.axesHandle;

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

plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 01.0], 'MarkerSize', 10);
axis('square')
set(gca, 'XLim', [otfRange/100 otfRange], 'XScale', 'log', 'YLim', [0 1.05]);
set(gca, 'XTick', sfTicks, 'YTick', 0:0.1:1.0);
set(gca, 'YTickLabel', 0:0.1:1);
grid('on'); box('on');
set(gca, 'FontSize', fontSize);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'bold'); 
ylabel('modulation');
title(sprintf('MTF (%2.0f nm)', targetWavelength));
end

