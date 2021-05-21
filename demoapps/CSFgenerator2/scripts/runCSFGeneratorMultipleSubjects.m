function runCSFGeneratorMultipleSubjects

    compute = true;
    if (compute)
        computeCSF()
    else
        visualizeCSF()
    end
    
end

function visualizeCSF()
    load('subjectsData.mat', 'sfSupport', ...
        'constantSizeWatsonPyradidOfVisbility', ...
        'constantCyclesWatsonPyradidOfVisbility', ...
        'csfDataConstantCycles', ...
        'csfDataConstantSize');
    
    % Display the computed CSF
    hFig = figure(1); clf;
    set(hFig, 'Position', [100 100 600 700], 'Color', [1 1 1]);
    
    % Superimpose Watson's pyradid of visibility data for constant size and
    % constant cycles
    plot(sfSupport, constantSizeWatsonPyradidOfVisbility, ...
        'b--', 'LineWidth', 1.5);
    hold on;
    plot(sfSupport, constantCyclesWatsonPyradidOfVisbility, ...
        'r--', 'LineWidth', 1.5);
    
    for subjectID = 1:10
        % constant size CSF
        csfDataCS = csfDataConstantSize{subjectID};

        % constant cycles CSF
        csfDataCC = csfDataConstantCycles{subjectID};
        
        plot(csfDataCS.spatialFrequencySupport, csfDataCS.sensitivity, ...
            'bo-', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.5 1.0], 'LineWidth', 1.5);
        plot(csfDataCC.spatialFrequencySupport, csfDataCC.sensitivity, ...
            'ro-', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5],'LineWidth', 1.5);
    
    end
    
    
    set(gca, 'XScale', 'linear', 'YScale', 'log', ...
        'XLim', [2 50], 'YLim', [1 300], 'FontSize', 16);
    set(gca, 'XTick', 0:5:40, ...
        'YTick', [1 3 10 30 100 300 1000 3000]);
    xlabel('spatial frequency (c/deg)');
    ylabel('sensitivity');
    grid on;
    legend({'Watson''s PofV constant size',  'Watson''s PoV constant cycles', 'ISETBio constant size', 'ISETBio constant cycles'});
    
    
end

function displayData()

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
    plot(sfSupport, constantSizeWatsonPyradidOfVisbility, ...
        'b--', 'LineWidth', 1.5);
    plot(sfSupport, constantCyclesWatsonPyradidOfVisbility, ...
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



function computeCSF()
    % Load params exported by the app
    load('CSFGeneratorAppExport.mat', 'params');
        
    % Display the entire params struct
    CSFGeneratorApp.render.paramsStructTree(params);
    
    for subjectID = 1:10
    
        params.optics.subjectID = subjectID;
        % Modify params here to run a constant-size stimulus
        params.csfParams.constantParameter = 'constant size';
        params.csfParams.spatialFrequencyMin = 3;
        params.csfParams.spatialFrequencyMax = 40;
        params.csfParams.spatialFrequencySamples = 7;

        params.stimParams.sizeDegs = 2.0;
        params.stimParams.resolutionPixels = 300;
        
        params.stimParams.meanLuminanceCdM2 = 30;
        params.coneMosaicParams.integrationTime = 125/1000;

        params.coneMosaicParams.sizeDegs = [2.0 2.0];
        params.coneMosaicParams.eccentricityDegs = [0 0];

        params.psychometricFunctionParams.testTrials = 300;
        params.debugParams.visualizeConeMosaicActivationComponents = ~true;

        % Build the non-gui part of the app
        app = params;
    
        % Generate components based on the params
        if (subjectID == 1)
            app.components.coneMosaic = CSFGeneratorApp.generate.coneMosaic(params, []);
            theConeMosaic = app.components.coneMosaic;
        else
            app.components.coneMosaic = theConeMosaic;
        end
        % Generate optics for the subject
        app.components.optics = CSFGeneratorApp.generate.optics(app, []);

        % Compute the CSF
        csfDataConstantSize{subjectID} = CSFGeneratorApp.compute.contrastSensitivityFunction(app);

        % Modify params to run a constant cycles stimulus
        app.csfParams.constantParameter = 'constant cycles';
        app.csfParams.numberOfConstantCycles = 0.56; % 0.56 cycles/sigma -> 1 octave bandwidth

        % Compute the CSF
        csfDataConstantCycles{subjectID} = CSFGeneratorApp.compute.contrastSensitivityFunction(app);
    end

    % Generate Watson's Pyramid of Visibility curve for constant size
    [sfSupport, constantSizeWatsonPyramidOfVisibility] = CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
        app.csfParams.spatialFrequencyMin,  app.csfParams.spatialFrequencyMax, ...
        app.stimParams.meanLuminanceCdM2, 'constant size');
    
    % Generate Watson's Pyramid of Visibility curve for constant cycles
    [sfSupport, constantCyclesWatsonPyramidOfVisibility] = CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
        app.csfParams.spatialFrequencyMin,  app.csfParams.spatialFrequencyMax, ...
        app.stimParams.meanLuminanceCdM2, 'constant cycles');
    
    
    save('subjectsData.mat', 'sfSupport', ...
        'constantSizeWatsonPyradidOfVisbility', ...
        'constantCyclesWatsonPyradidOfVisbility', ...
        'csfDataConstantCycles', ...
        'csfDataConstantSize');
    
end



