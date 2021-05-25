function [csfData, psychometricFunctionData, stopRun] = contrastSensitivityFunction(app, varargin)
    
    if (isa(app, 'ISETBioCSFGenerator'))
        appCall = true;
    else
        appCall = false;
    end

    % Form output struct
    csfData = struct(....
        'spatialFrequencySupport', logspace(...
                                        log10(app.csfParams.spatialFrequencyMin), ...
                                        log10(app.csfParams.spatialFrequencyMax), ...
                                        app.csfParams.spatialFrequencySamples), ...
        'sensitivity', zeros(1, app.csfParams.spatialFrequencySamples));
    

    % Instantiate a neural response engine using the
    % @nrePhotopigmentExcitationsCmosaicSingleShot compute function
    theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaicSingleShot);
    theNeuralEngine.customNeuralPipeline(struct(...
          'coneMosaic', app.components.coneMosaic, ...
          'optics', app.components.optics));
      
    % Instantiate a response classifier engine using the @rcePoissonTAFC
    theClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
      
    % Classifier engine params
    classifierEngineParams = struct(...
         'nTrain', 1, ...                                             % 1 training response
         'trainFlag', 'none', ...                                     % no noise for the single training response
         'nTest', app.psychometricFunctionParams.testTrials, ...      % test trials
         'testFlag', 'random'  ...                                    % random noise for the test responses
         );
        
    % Parameters for threshold estimation/quest engine
    thresholdParams = struct(...
        'logThreshLimitLow', -app.psychometricFunctionParams.log10ContrastMin, ...
        'logThreshLimitHigh', -app.psychometricFunctionParams.log10ContrastMax, ...
        'logThreshLimitDelta', app.psychometricFunctionParams.log10ContrastDelta, ...
        'slopeRangeLow', app.psychometricFunctionParams.slopeMin, ...
        'slopeRangeHigh', app.psychometricFunctionParams.slopeMax, ...
        'slopeDelta', app.psychometricFunctionParams.slopeDelta);
            
    % Quest engine parameters
    trialsNum = app.psychometricFunctionParams.testTrials * app.psychometricFunctionParams.contrastLevels;
    questEngineParams = struct(...
        'minTrial', trialsNum, 'maxTrial', trialsNum, ...
        'numEstimator', 1, 'stopCriterion', 0.05);
            
    % Generate stimParams for gratingSceneEngine
    [fixedStimParamsStruct, coneMosaicIsTooSmall] = CSFGeneratorApp.generate.stimParamsStructForGratingSceneEngine(app, app.stimParams.spatialFrequencyCPD);
    
    if (coneMosaicIsTooSmall)
        message = sprintf('Cone mosaic (%2.2f x %2.2f degs) does not extend over the stimulus (%2.2f degs). Increase size of cone mosaic or decrease size of stimulus.\n', ...
            app.coneMosaicParams.sizeDegs(1), app.coneMosaicParams.sizeDegs(2), CSFGeneratorApp.compute.stimulusSizeAtLowestSpatialFrequency(app));
        if (appCall)
            % Some feecback
            app.statusMessages('computational observer')  = struct(...
                'text', message,...
                'fontColor', app.colors('problem message foreground'), ...
                'backgroundColor', app.colors('problem message background'), ...
                'fontWeight', 'normal');

            % Render the status on the status field of tab B
            CSFGeneratorApp.render.statusField(app,'B', 'computational observer'); 
            return;
        else
            error(message);
        end
    end
    
    app.resetButtonPressed = false;
    stopRun = false;
    
    tic
    % Go through each spatial frequency and estimate the psychometric function
    for iSF = 1:length(csfData.spatialFrequencySupport) 
        
        if (stopRun)
            continue;
        end
        
        % Spatial frequency examined
        [fixedStimParamsStruct, coneMosaicIsTooSmall] = CSFGeneratorApp.generate.stimParamsStructForGratingSceneEngine(app, csfData.spatialFrequencySupport(iSF));
        
        % Contrast will be varied by mQuest
        fixedStimParamsStruct.contrast = [];
        
        % Throw an error if stimulus is out of Gamut
        fixedStimParamsStruct.warningInsteadOfErrorOnOutOfGamut = false;

        % Some feedback
        message = sprintf('Computing psychometric function for %2.1f c/deg\n', fixedStimParamsStruct.sf);
        if (appCall)
            if (app.resetButtonPressed)
                message = 'Run aborted';
                stopRun = true;
            end
            % Some feecback
            app.statusMessages('computational observer')  = struct(...
                'text', message, ...
                'fontColor', app.colors('good message foreground'), ...
                'backgroundColor', app.colors('good message background'), ...
                'fontWeight', 'normal');
            % Render the status on the status field of tab B
            CSFGeneratorApp.render.statusField(app,'B', 'computational observer');
        else
            fprintf('%s\n', message);
        end
        

        % Generate the stimulus scene engine
        [~, ~, ~, gratingSceneEngine] = CSFGeneratorApp.generate.gratingSceneEngine(fixedStimParamsStruct, []);
        
        % Compute the threshold for this grating  using the imported neural
        % and classifier engines and params
        [logThresholdData, questOBJ] = computeThresholdTAFC(...
            gratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
            classifierEngineParams, ...
            thresholdParams, ...
            questEngineParams, ...
            'visualizeAllComponents', app.debugParams.visualizeConeMosaicActivationComponents, ...
            'beVerbose', false);
        
        % Call the thresholdMLE() method of the questOBJ to update
        % the current psychometric function data
        [~, ~, psychometricFunctionData{iSF}] = questOBJ.thresholdMLE(...
                    'showPlot',  false, 'newFigure', false, ...
                    'returnData', true);
        
        % Save data
        thePsychometricFunction = psychometricFunctionData{iSF};
        thresholdPCorrect = 0.81;
        if (max(thePsychometricFunction.pCorrect) > thresholdPCorrect)
            csfData.sensitivity(iSF) = 10^(-logThresholdData); 
        else
            csfData.sensitivity(iSF) = nan;
        end
        
        if (appCall)
            % Plot the estimated psychometric function for this spatial frequency
            CSFGeneratorApp.render.psychometricFunctionView(app, 'update', ...
                'withData', psychometricFunctionData);
            
            % Display the computed CSF
            CSFGeneratorApp.render.csfView(app, 'update', 'withData', csfData);

        end
    end % iSF
    
    app.resetButtonPressed = false;
    
    if (appCall)
        % Display the computed CSF
        CSFGeneratorApp.render.csfView(app, 'update', 'withData', csfData);

        if (stopRun)
            message = 'Run aborted by user';
        else
            message = sprintf('Contrast sensitivity function computation completed in %2.1f minutes.', toc/60);
        end
        % Some feecback
        app.statusMessages('computational observer')  = struct(...
                'text', message,...
                'fontColor', app.colors('good message foreground'), ...
                'backgroundColor', app.colors('good message background'), ...
                'fontWeight', 'normal');

        % Render the status on the status field of tab B
        CSFGeneratorApp.render.statusField(app,'B', 'computational observer'); 
    end
    
    % Save the data
    app.products.psychometricFunctionData = psychometricFunctionData;
    app.products.csfData = csfData;
end