function runCSFdemo()

    % Load pipeline data exported from the CSF app
    load('CSFappExport.mat');
    
    % Define the test spatial frequencies
    examinedSpatialFrequencies = [2 4 8 16 32 64];
    
    % Compute the CSF for the test spatial frequencies using the data exported from the app
    csf = CSFgenerator.computeCSF(examinedSpatialFrequencies, exportedData.stimData, exportedData.pipeline);
    
    % Load another pipeline also exported from the CSF app.
    % In this pipeline optics are corrected for the central defocus
    load('CSFappExportSubtractCentralRefraction.mat');
    csfRefracted = CSFgenerator.computeCSF(examinedSpatialFrequencies, exportedData.stimData, exportedData.pipeline);
    
    % Plot data
    hFig = figure(1);
    set(hFig, 'Position', [100 100 500 800]);
    plot(examinedSpatialFrequencies, csf, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]); hold on;
    plot(examinedSpatialFrequencies, csfRefracted, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    legend({'uncorrected foveal defocus', 'zero foveal defocus'});
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [1 100], 'YLim', [1 1000]);
    set(gca, 'FontSize', 14);
    xlabel('spatial frequency (c/deg)');
    ylabel('sensitivity');
    
end
