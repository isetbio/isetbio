function runCSFGenerator

    % Load params exported by the app
    load('CSFGeneratorAppExport.mat', 'params');
        
    % Display the entire params struct
    CSFGeneratorApp.render.paramsStructTree(params);
    pause;
    
    % Modify params here to run a constant-size stimulus
    params.csfParams.constantParameter = 'constant size';
    params.csfParams.spatialFrequencyMin = 3;
    params.csfParams.spatialFrequencyMax = 40;
    params.csfParams.spatialFrequencySamples = 6;
    
    params.stimParams.sizeDegs = 2.0;
    params.stimParams.resolutionPixels = 300;
    
    params.stimParams.meanLuminanceCdM2 = 30;
    
    % 125 msec integration time: single-shot, no eye movements
    params.coneMosaicParams.integrationTime = 125/1000;
    
    % Generate a 2x2 cone mosaic
    params.coneMosaicParams.sizeDegs = [2.0 2.0];
    
    % Centered at the origin
    params.coneMosaicParams.eccentricityDegs = [0 0];
    
    params.psychometricFunctionParams.testTrials = 512;
    params.debugParams.visualizeConeMosaicActivationComponents = ~true;
    
    % Build the non-gui part of the app
    app = params;
    
    % Generate components based on the params
    app.components.coneMosaic = CSFGeneratorApp.generate.coneMosaic(params, []);
    app.components.optics = CSFGeneratorApp.generate.optics(app, []);
    
    % Compute the CSF for constant size stimuli
    app.csfParams.constantParameter = 'constant size';
    csfDataConstantSize = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
     
    % % Compute the CSF for constant cycles stimuli
    app.csfParams.constantParameter = 'constant cycles';
    app.csfParams.numberOfConstantCycles = 0.56; % 0.56 cycles/sigma -> 1 octave bandwidth
    csfDataConstantCycles = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
    
    % Generate Watson's Pyramid of Visibility curve for constant size
    [sfSupport, constantSizeWatsonPyramidOfVisibility] = CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
        params.csfParams.spatialFrequencyMin,  params.csfParams.spatialFrequencyMax, ...
        app.stimParams.meanLuminanceCdM2, 'constant size');
    
    % Generate Watson's Pyramid of Visibility curve for constant cycles
    [sfSupport, constantCyclesWatsonPyramidOfVisibility] = CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
        params.csfParams.spatialFrequencyMin,  params.csfParams.spatialFrequencyMax, ...
        app.stimParams.meanLuminanceCdM2, 'constant cycles');
    
    % Display the computed CSF
    hFig = figure(1); clf;
    set(hFig, 'Position', [100 100 600 700], 'Color', [1 1 1]);
    
    plot(csfDataConstantSize.spatialFrequencySupport, csfDataConstantSize.sensitivity, ...
        'bo-', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.5 1.0], 'LineWidth', 1.5);
    hold on
    plot(csfDataConstantCycles.spatialFrequencySupport, csfDataConstantCycles.sensitivity, ...
        'ro-', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5],'LineWidth', 1.5);
   
    % Superimpose Watson's pyradid of visibility data for constant size and
    % constant cycles
    plot(sfSupport, constantSizeWatsonPyramidOfVisibility, ...
        'b--', 'LineWidth', 1.5);
    plot(sfSupport, constantCyclesWatsonPyramidOfVisibility, ...
        'r--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'linear', 'YScale', 'log', ...
        'XLim', [2 70], 'YLim', [1 300], 'FontSize', 16);
    set(gca, 'XTick', [2 3 5 10 20 30 50], ...
        'YTick', [1 3 10 30 100 300 1000 3000]);
    xlabel('spatial frequency (c/deg)');
    ylabel('sensitivity');
    grid on;
    legend({'ISETBio constant size', 'ISETBio constant cycles', 'Watson''s PofV constant size', 'Watson''s PoV constant cycles'});
    
end


