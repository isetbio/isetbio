function [theNeuralEngine, theClassifierEngine] = generateNeuralAndClassifierEngines(pipelineParams)

    % Generate the neural engine
    theNeuralEngine = neuralResponseEngine(pipelineParams.neuralComputeFunction);
    theNeuralEngine.customNeuralPipeline(struct(...
          'coneMosaic', pipelineParams.coneMosaic, ...
          'optics', pipelineParams.optics));
    
    % Generate the classifier engine
    theClassifierEngine = responseClassifierEngine(pipelineParams.classifierComputeFunction);
                    
end

