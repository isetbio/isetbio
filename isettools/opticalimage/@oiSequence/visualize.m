% Visualize the sequence
function visualize(obj)

colsNum = round(1.3*sqrt(obj.length));
rowsNum = ceil(obj.length/colsNum);
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum, ...
           'colsNum', colsNum+1, ...
           'heightMargin',   0.07, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);

XYZmax = 0;
for oiIndex = 1:obj.length
     currentOI = obj.frameAtIndex(oiIndex);
     XYZ = oiGet(currentOI, 'xyz');
     if (max(XYZ(:)) > XYZmax)
         XYZmax = max(XYZ(:));
     end
end
% Do not exceed XYZ values of 0.5 (for correct rendering)
XYZmax = 2*XYZmax;

hFig = figure();
set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 730]);

for oiIndex = 1:obj.length
    if (oiIndex == 1)
        % Plot the modulation function
        subplot('Position', subplotPosVectors(1,1).v);
        stairs(1:obj.length, obj.modulationFunction, 'r', 'LineWidth', 1.5);
        set(gca, 'XLim', [1 obj.length], 'FontSize', 12);
        title(sprintf('composition: ''%s''', obj.composition));
        xlabel('frame index');
        ylabel('modulation');
    end
    
    % Ask theOIsequence to return the oiIndex-th frame
    currentOI = obj.frameAtIndex(oiIndex);
    support = oiGet(currentOI, 'spatial support', 'microns');
    [~, meanIlluminance] = oiCalculateIlluminance(currentOI);
    xaxis = support(1,:,1);
    yaxis = support(:,1,2);
    row = 1+floor((oiIndex)/(colsNum+1));
    col = 1+mod((oiIndex),(colsNum+1));
    
    subplot('Position', subplotPosVectors(row,col).v);
    rgbImage = xyz2srgb(oiGet(currentOI, 'xyz')/XYZmax);
    imagesc(xaxis, yaxis, rgbImage, [0 1]);
    axis 'image'
    if (col == 1) && (row == rowsNum)
       xticks = [xaxis(1) 0 xaxis(end)];
       yticks = [yaxis(1) 0 yaxis(end)];
       set(gca, 'XTick', xticks, 'YTick', yticks, 'XTickLabel', sprintf('%2.0f\n', xticks), 'YTickLabel', sprintf('%2.0f\n', yticks));
       ylabel('microns');
    else
       set(gca, 'XTick', [], 'YTick', [])
       xlabel(sprintf('frame %d (%2.1fms)', oiIndex, 1000*obj.oiTimeAxis(oiIndex)));
    end
    title(sprintf('mean illum: %2.1f', meanIlluminance));
    set(gca, 'FontSize', 12);
end
 
end

