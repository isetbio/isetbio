function runCSFGenerator

    % Load params exported by the app
    load('CSFGeneratorAppExport.mat', 'params');
        
    % Display the entire params struct
    CSFGeneratorApp.render.paramsStructTree(params);
    pause
    
    % Modify params here
    params.stimParams.sizeDegs = 1;
    params.stimParams.mosaicCenteredPosition = true;
    params.csfParams.constantParameter = 'constant size';
    params.coneMosaicParams.sizeDegs = [1.2 1.2];
    params.coneMosaicParams.eccentricityDegs = [0 0];
    params.debugParams.visualizeConeMosaicActivationComponents = ~true;
    
   % params.csfParams.constantParameter = 'constant cycles';
   % params.csfParams.numberOfConstantCycles = 0.56;
    ...
    
    % Build the non-gui part of the app
    app = params;
    
    % Generate components based on the params
    app.components.coneMosaic = CSFGeneratorApp.generate.coneMosaic(params, []);
    app.components.optics = CSFGeneratorApp.generate.optics(app, []);
    
    % Compute the CSF
    csfData = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
     
    % Display the computed CSF
    figure(1); hold on
    plot(csfData.spatialFrequencySupport, log10(csfData.sensitivity), 'ro-', 'MarkerSize', 16, 'LineWidth', 1.0);
    set(gca, 'XScale', 'linear', 'YScale', 'linear', 'XLim', [0 40], 'YLim', [0 4], 'FontSize', 16);
end


