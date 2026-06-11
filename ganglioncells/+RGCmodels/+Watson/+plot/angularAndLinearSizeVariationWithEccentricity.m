function angularAndLinearSizeVariationWithEccentricity()
    
    % We will examine a stimulus whose linear size is 50 microns on the retina
    constantSizeRetinalMicrons = 100; 
    % And a stimulus whose angular size is 1 visual degree
    constantSizeVisualDegs = 1.0;
    
    % Eccentricity range to examine
    eccMinDegs = 0.0;
    eccMaxDegs = 100;
    eccSamplesNum = 200;
    
    eccDegs = linspace(eccMinDegs, eccMaxDegs, eccSamplesNum);
    eccMicrons = RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs)*1e3;
    
    % Convert retinal size in microns to visual degrees
    sizeDegs = RGCmodels.Watson.convert.sizeRetinalMicronsToSizeVisualDegs(constantSizeRetinalMicrons, eccMicrons);
    
    % Convert visual size in degrees to retinal microns
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(constantSizeVisualDegs, eccDegs);
    
    hFig = figure(99); clf;
    set(hFig, 'Position', [10 10 1000 500], 'Color', [1 1 1]);
    
    subplot(1,2,1);
    plot(eccDegs, sizeDegs, 'r-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0.3 100], 'XScale', 'log', 'XTick', [0.3 1 3 10 30 100], ...
        'YLim', [0.2 0.8], 'FontSize', 14);
    axis 'square'; grid on
    xlabel('eccentricity (degs)');
    ylabel('angular size (visual degs)');
    title(sprintf('angular size of a %2.0f micron stimulus', constantSizeRetinalMicrons));
    
    subplot(1,2,2);
    plot(eccDegs, sizeMicrons, 'r-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0.3 100], 'XScale', 'log', 'XTick', [0.3 1 3 10 30 100], ...
        'YLim', [200 300], 'FontSize', 14);
    axis 'square'; grid on
    xlabel('eccentricity (degs)');
    ylabel('linear size (retinal microns)');
    title(sprintf('retinal size of a %2.0f deg stimulus', constantSizeVisualDegs));
    
    
    
end

