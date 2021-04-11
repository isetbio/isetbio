function csf = computeCSF(sfSupport, stimParams, pipelineParams)

    % Generate neural and classifier engines
    [theNeuralEngine, theClassifierEngine] = CSFgenerator.generateNeuralAndClassifierEngines(pipelineParams);
    
    csf = zeros(1, numel(sfSupport));
    for idx = 1:length(sfSupport) 
        % Some feedback
        fprintf('Computing sensitivity to %2.1f c/deg\n', sfSupport(idx));
        
        % Call the external CSFgenerator.computeStimulusSceneEngine method to
        % generate the gratingSceneEngine
        sParams = stimParams;
        sParams.sf = sfSupport(idx);
        [~, ~, ~, gratingSceneEngine] = CSFgenerator.computeStimulusSceneEngine(sParams, []);
              
        % Compute the threshold for this grating  using the imported neural
        % and classifier engines and params
        logThresholdData = computeThresholdTAFC(...
            gratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
            pipelineParams.classifierEngineParams, ...
            pipelineParams.thresholdParams, ...
            pipelineParams.questEngineParams, ...
            'beVerbose', false);
                
        csf(idx) = 10^(-logThresholdData); 
    end % idx
    
end

