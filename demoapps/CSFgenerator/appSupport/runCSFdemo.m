function runCSFdemo()

    % Load stimulus, pipeline and CSF data exported from the CSF app
    load('CSFappExport.mat', 'exportedData');
   
    % Extract the csf data computed by the app
    [sfOriginal, csfOriginal] = extractImportedCSFdata(exportedData.csfData);

    % Extract the stimulus data from the app
    stimData = exportedData.stimData;

    % Replicate the csf computed by the app (for comparison)
    csfRepeat = CSFgenerator.computeCSF(sfOriginal, stimData, exportedData.pipeline, ...
        'visualizeAllComponents', false);
    
    % Increase the stimulus size to 1.5 degs
    stimData.fovDegs = 1.5;
    stimData.spatialEnvelopeRadiusDegs = stimData.fovDegs/2;
    
    % Recompute CSF at these spatial frequencies
    examinedSpatialFrequencies = [2 4 8 16 32 64];
    % Go !
    csf = CSFgenerator.computeCSF(examinedSpatialFrequencies, stimData, exportedData.pipeline, ...
        'visualizeAllComponents', false);
    
    % Plot data
    hFig = figure(1);
    set(hFig, 'Position', [100 100 500 800]);
    plot(sfOriginal, csfOriginal, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]); hold on;
    plot(sfOriginal, csfRepeat, 'ks-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    plot(examinedSpatialFrequencies, csf, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    legend({'from app', 'repeated', 'recomputed using a 1.5 deg stimulus'});
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.5 100], 'YLim', [1 1000]);
    set(gca, 'XTick', [1 3 10 30 60 100], 'YTick', [1 3 10 30 100 300 1000 3000]);
    grid 'on'
    set(gca, 'FontSize', 14);
    xlabel('spatial frequency (c/deg)');
    ylabel('sensitivity');
    
end

function [sf, csf] = extractImportedCSFdata(csfData)

    sf = zeros(1, length(csfData.csf));
    csf = zeros(1, length(csfData.csf));
    
    for i = 1:length(csfData.csf)
        sf(i) = csfData.csf{i}.spatialFrequency;
        csf(i) = csfData.csf{i}.sensitivity;
    end
    
end
