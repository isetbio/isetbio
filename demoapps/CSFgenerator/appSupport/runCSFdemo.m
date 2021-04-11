function runCSFdemo()

    % Load pipeline data exported from the CSF app
    load('CSFappExport.mat');
    
    % Define the test spatial frequencies
    examinedSpatialFrequencies = [2 4 8 16 32 64];
    
    % Compute the CSF for the test spatial frequencies using the data exported from the app
    csf = CSFgenerator.computeCSF(examinedSpatialFrequencies, exportedData.stimData, exportedData.pipeline);
    
    load('CSFappExportSubtractCentralRefraction.mat');
    csfRefracted = CSFgenerator.computeCSF(examinedSpatialFrequencies, exportedData.stimData, exportedData.pipeline);
    
    % Plot data
    figure();
    plot(examinedSpatialFrequencies, csf, 'ko-'); hold on;
    plot(examinedSpatialFrequencies, csfRefracted, 'ro-');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [1 100], 'YLim', [1 1000]);
end
