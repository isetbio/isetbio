function runCSFGenerator

    % Load params exported by the app
    load('CSFGeneratorAppExport.mat', 'params');
        
    % Display the entire params struct
    CSFGeneratorApp.render.paramsStructTree(params);
    pause
    
    % Modify params here
    ...
    
    % Build the non-gui part of the app
    app = params;
    
    % Generate components based on the params
    app.components.coneMosaic = CSFGeneratorApp.generate.coneMosaic(params, []);
    app.components.optics = CSFGeneratorApp.generate.optics(app, []);
    
    % Compute the CSF
    csfData = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
     
    % Display the computed CSF
    figure(1);
    plot(csfData.spatialFrequencySupport, csfData.sensitivity, 'ko-');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [1 100], 'FontSize', 16);
end


