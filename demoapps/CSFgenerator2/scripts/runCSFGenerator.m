function runCSFGenerator

    % Load params exported by the app
    load('CSFGeneratorAppExport.mat', 'params');
        
    % Display the entire params struct
    CSFGeneratorApp.render.paramsStructTree(params);
    
    % Modify params here to run a constant-size stimulus
    params.csfParams.constantParameter = 'constant size';
    params.csfParams.spatialFrequencyMin = 3;
    params.csfParams.spatialFrequencyMax = 60;
    params.csfParams.spatialFrequencySamples = 7;
    
    params.stimParams.sizeDegs = 2.5;
    params.stimParams.resolutionPixels = 400;
    
    params.stimParams.meanLuminanceCdM2 = 30;
    params.coneMosaicParams.integrationTime = 125/1000;
    
    params.coneMosaicParams.sizeDegs = [2.5 2.5];
    params.coneMosaicParams.eccentricityDegs = [0 0];
    
    params.psychometricFunctionParams.testTrials = 512;
    params.debugParams.visualizeConeMosaicActivationComponents = ~true;
    
    % Build the non-gui part of the app
    app = params;
    
    % Generate components based on the params
    app.components.coneMosaic = CSFGeneratorApp.generate.coneMosaic(params, []);
    app.components.optics = CSFGeneratorApp.generate.optics(app, []);
    
    % Compute the CSF
    csfDataConstantSize = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
     
    % Modify params to run a constant cycles stimulus
    app.csfParams.constantParameter = 'constant cycles';
    app.csfParams.numberOfConstantCycles = 0.56; % 0.56 cycles/sigma -> 1 octave bandwidth

    % Compute the CSF
    csfDataConstantCycles = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
    
    % Generate Watson's Pyramid of Visibility curve for constant size
    [sfSupport, constantSizeWatsonPyradidOfVisbility] = WatsonData(app, 'CDG');
    
    % Generate Watson's Pyramid of Visibility curve for constant cycles
    [sfSupport, constantCyclesWatsonPyradidOfVisbility] = WatsonData(app, 'CCG');
    
    % Display the computed CSF
    figure(1); clf;
    plot(csfDataConstantSize.spatialFrequencySupport, csfDataConstantSize.sensitivity, ...
        'bo-', 'MarkerSize', 16, 'LineWidth', 1.5);
    hold on
    plot(csfDataConstantCycles.spatialFrequencySupport, csfDataConstantCycles.sensitivity, ...
        'ro-', 'MarkerSize', 16, 'LineWidth', 1.5);
   
    % Superimpose Watson's pyradid of visibility data for constant size and
    % constant cycles
    plot(sfSupport, constantSizeWatsonPyradidOfVisbility, ...
        'b--', 'MarkerSize', 16, 'LineWidth', 1.5);
    plot(sfSupport, constantCyclesWatsonPyradidOfVisbility, ...
        'r--', 'MarkerSize', 16, 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0 40], 'YLim', [0 2000], 'FontSize', 16);
    set(gca, 'XTick', [3 10 30], 'YTick', [1 3 10 30 100 300 1000 3000]);
    grid on;
    legend({'ISETBio constant size', 'ISETBio constant cycles', 'Watson constant size', 'Watson constant cycles'});
    
end


function [sfSupport, S] = WatsonData(app, stimSize)
    temporalFrequency = 0.0; cW = 0;
    logLuminanceNits = log10(app.stimParams.meanLuminanceCdM2);
    % Table 1, of Watson 2018, "The Field of View, the Field of Resolution and the 
    %       Field of Contrast Sensitivity" (Luminance)
    switch (stimSize)
        case 'CDG'
            c0 = 1.739; cF = -0.060; cL = 0.391;
        case 'CCG'
            c0 = 1.380; cF = -0.091; cL = 0.391;
    end
    
    sfSupport = 3:1:60;
    logS = c0 + cW*temporalFrequency + cF*sfSupport + cL*logLuminanceNits;
    S = 10.^logS;
end

