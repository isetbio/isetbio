function hFig = visualizePSFsubStack(theOI, sampledWavelengths, psfRangeArcMin)
  
    if (numel(sampledWavelengths)>10)
        rowsNum = 3;
    elseif (numel(sampledWavelengths)>5)
        rowsNum = 2;
    else
        rowsNum = 1;
    end

    colsNum = ceil(numel(sampledWavelengths)/rowsNum);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', rowsNum, ...
        'colsNum', colsNum, ...
        'heightMargin',  0.04, ...
        'widthMargin',    0.01, ...
        'leftMargin',     0.04, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.07, ...
        'topMargin',      0.02);
    
    wavelengthSupport = oiGet(theOI, 'wave');

    hFig = figure(); clf;
    set(hFig,  'Color', [1 1 1], 'Position', [10 10 1600 900]);

    psfColorMap = brewermap(1024, 'greys');

    for k = 1:numel(sampledWavelengths)
   
        [~,wIndex] = min(abs(wavelengthSupport-sampledWavelengths(k)));
        targetWavelength = wavelengthSupport(wIndex);

        row = floor((k-1)/colsNum)+1;
        col = mod(k-1,colsNum)+1;
        
        noXLabel = true;
        noYLabel = true;
        if (row == rowsNum)
            noXLabel = false;
        end
        if (col == 1)
            noYLabel = false;
        end
        ax = subplot('Position', subplotPosVectors(row,col).v);
        visualizePSF(theOI, targetWavelength, psfRangeArcMin, ...
            'contourLevels', 0.1:0.1:0.9, ...
            'axesHandle', ax, ...
            'psfColorMap', psfColorMap, ...
            'noXLabel', noXLabel, ...
            'noYLabel', noYLabel, ...
            'figureTitle', sprintf('%2.0f nm', targetWavelength), ...
            'fontSize', 20);
    end
    
end