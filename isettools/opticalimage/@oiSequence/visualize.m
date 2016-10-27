% Visualize the sequence
function visualize(obj)

colsNum = round(sqrt(obj.length));
rowsNum = floor(obj.length/colsNum);
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum, ...
           'colsNum', colsNum+1, ...
           'heightMargin',   0.025, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
       
for oiIndex = 1:obj.length
    if (oiIndex == 1)
        % Plot the modulation function
        subplot('Position', subplotPosVectors(1,1).v);
        plot(1:obj.length, obj.modulationFunction, 'rs-', 'LineWidth', 1.5);
        set(gca, 'XLim', [1 obj.length], 'FontSize', 12);
        title(sprintf('composition\n''%s''', obj.composition));
        xlabel('frame index');
        ylabel('modulation');
    end
    
    % Ask theOIsequence to return the oiIndex-th frame
    currentOI = obj.frameAtIndex(oiIndex);
    support = oiGet(currentOI, 'spatial support', 'microns');
    xaxis = support(1,:,1);
    yaxis = support(:,1,2);
    row = 1+floor((oiIndex)/(colsNum+1));
    col = 1+mod((oiIndex),(colsNum+1));
    
    subplot('Position', subplotPosVectors(row,col).v);
    rgbImage = xyz2srgb(oiGet(currentOI, 'xyz'));
    imagesc(xaxis, yaxis, rgbImage, [0 1]);
    if (col == 1) && (row == rowsNum)
       xticks = [xaxis(1) 0 xaxis(end)];
       yticks = [yaxis(1) 0 yaxis(end)];
       set(gca, 'XTick', xticks, 'YTick', yticks, 'XTickLabel', sprintf('%2.0f\n', xticks), 'YTickLabel', sprintf('%2.0f\n', yticks));
    else
       set(gca, 'XTick', [], 'YTick', [])
       xlabel(sprintf('frame %d', oiIndex));
    end

    axis 'image'
    set(gca, 'FontSize', 12);
end
 
end

