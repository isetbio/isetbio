function visualizeMidgetMosaicsAtDifferentEccentricities

    eccSizeDegsExamined = [...
         0  0.5; ...
        -1  0.5; ...
        -2  0.5; ...
        -3  0.6; ...
        -4  0.8; ...
        -6  1.0; ...
        -8  1.2; ...
        -10 1.4; ...
        -12 1.6; ...
        -14 1.8; ...
        -16 2.0; ...
        -20 2.2; ...
        -25 0.5 ...
       ];
    
    mRGCMosaicIndex = size(eccSizeDegsExamined,1);
    sizeDegs = [0.7 0.7];

    
    for k = 12:12
        eccIndex = mRGCMosaicIndex-k;

        theMidgetRGCmosaic = midgetRGCMosaic(...
           'sourceLatticeSizeDegs', 60, ...
           'eccentricityDegs', [eccSizeDegsExamined(eccIndex,1) 0], ...
           'sizeDegs', sizeDegs ...
        );
    
        hFig = figure(1); clf;
        set(hFig, 'Position', [100 100 700 700], 'Color', [1 1 1]);
        ax = subplot('Position', [0.12 0.08 0.85 0.91]);
        theMidgetRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'maxVisualizedRFs', 22, ...
            'yLims', [-0.15 0.15], ...
            'fontSize', 20);
    
        NicePlot.exportFigToPDF(sprintf('Mosaic_%2.1fdegs.pdf',eccSizeDegsExamined(eccIndex,1)), hFig, 300);
    end

end