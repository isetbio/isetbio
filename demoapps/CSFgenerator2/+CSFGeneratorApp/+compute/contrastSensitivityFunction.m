function contrastSensitivityFunction = contrastSensitivityFunction(app, varargin)
    
    p = inputParser;
    p.addParameter('appIsParams', false, @islogical);
    p.parse(varargin{:});
    
    % Form output struct
    contrastSensitivityFunction = struct(....
        'spatialFrequencySupport', logspace(...
                                        log10(app.csfParams.spatialFrequencyMin), ...
                                        log10(app.csfParams.spatialFrequencyMax), ...
                                        app.csfParams.spatialFrequencySamples), ...
        'sensitivity', zeros(1, app.csfParams.spatialFrequencySamples));
    


    theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaicSingleShot);
    theNeuralEngine.customNeuralPipeline(struct(...
          'coneMosaic', app.components.coneMosaic, ...
          'optics', app.components.optics));
      
    theClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
      
    % Classifier engine params
    classifierEngineParams = struct(...
         'nTrain', 1, ...                                             % 1 training response
         'trainFlag', 'none', ...                                     % no noise for the single training response
         'nTest', app.psychometricFunctionParams.testTrials, ...          % test trials
         'testFlag', 'random'  ...                                    % random noise for the test responses
         );
        
    % Parameters for threshold estimation/quest engine
    thresholdParams = struct(...
        'logThreshLimitLow', -app.psychometricFunctionParams.log10ContrastMin, ...
        'logThreshLimitHigh', -app.psychometricFunctionParams.log10ContrastMax, ...
        'logThreshLimitDelta', app.psychometricFunctionParams.log10ContrastDelta, ...
        'slopeRangeLow', app.psychometricFunctionParams.slopeMin, ...
        'slopeRangeHigh', app.psychometricFunctionParams.slopeMax, ...
        'slopeDelta', app.psychometricFunctionParams.slopeDelta)
            
    % Quest engine parameters
    trialsNum = app.psychometricFunctionParams.testTrials * app.psychometricFunctionParams.contrastLevels;
    questEngineParams = struct(...
        'minTrial', trialsNum, 'maxTrial', trialsNum, ...
        'numEstimator', 1, 'stopCriterion', 0.05)
            
    % Generate stimParams for gratingSceneEngine
    fixedStimParamsStruct = CSFGeneratorApp.generate.stimParamsStructForGratingSceneEngine(app.stimParams);
          
    visualizeAllComponents = true;
    for idx = 1:length(contrastSensitivityFunction.spatialFrequencySupport) 

        % Contrast will be varied by mQuest
        fixedStimParamsStruct.contrast = [];
        
        % Throw an error if stimulus is out of Gamut
        fixedStimParamsStruct.warningInsteadOfErrorOnOutOfGamut = false;
        
        % Spatial frequency examined
        fixedStimParamsStruct.sf = contrastSensitivityFunction.spatialFrequencySupport(idx);
        
        % Generate the stimulus scene engine
        [~, ~, ~, gratingSceneEngine] = CSFGeneratorApp.generate.gratingSceneEngine(fixedStimParamsStruct, []);
        
        % Some feedback
        fprintf('Computing psychometric function for %2.1f c/deg\n', fixedStimParamsStruct.sf);
        
        % Compute the threshold for this grating  using the imported neural
        % and classifier engines and params
        logThresholdData = computeThresholdTAFC(...
            gratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
            classifierEngineParams, ...
            thresholdParams, ...
            questEngineParams, ...
            'visualizeAllComponents', visualizeAllComponents, ...
            'beVerbose', false);
                
        contrastSensitivityFunction.sensitivity(idx) = 10^(-logThresholdData); 
    end % idx
    
end
        
            
