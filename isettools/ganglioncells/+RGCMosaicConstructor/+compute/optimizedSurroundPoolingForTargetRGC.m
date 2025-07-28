function  dataOut = optimizedSurroundPoolingForTargetRGC(theMRGCMosaic,  ...
  	theTargetRGCindex, targetVisualSTFparams, rfSurroundConnectivityParams, ...
  	stimParams, inputConeMosaicSTFresponses, varargin)

  % Parse input
  p = inputParser;
  p.addParameter('userSuppliedInitialValuesForModelVariables', [], @(x)(isempty(x)||(isstruct(x))));
  p.addParameter('visualizeFullAndMaximalExcursionSTF', false, @islogical);
  p.addParameter('visualizeGaussianFitToCenterSTF', false, @islogical);

  p.parse(varargin{:});
  visualizeFullAndMaximalExcursionSTF = p.Results.visualizeFullAndMaximalExcursionSTF;
  visualizeGaussianFitToCenterSTF = p.Results.visualizeGaussianFitToCenterSTF;
  userSuppliedInitialValuesForModelVariables = p.Results.userSuppliedInitialValuesForModelVariables;

 
  % Compute visual Rc (degs) by fitting the centerSTF with a Gaussian
  [visualRcDegs, theOptimalOrientation, anatomicalRcDegs] = ...
    RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.visualRcAtOptimalOrientation(...
      theMRGCMosaic, theTargetRGCindex, inputConeMosaicSTFresponses, stimParams, ...
      rfSurroundConnectivityParams.fitParams.GaussianToSTF.multiStartsNum, ...
      rfSurroundConnectivityParams.fitParams.GaussianToSTF.deltaThresholdToDeclareLocalMinInSTF, ...
      visualizeFullAndMaximalExcursionSTF, visualizeGaussianFitToCenterSTF, ...
      'fixedOptimalOrientation', rfSurroundConnectivityParams.fixedOptimalOrientation);
  

  % Assemble modelConstants struct
  modelConstants = struct(...
    'poolingModel', rfSurroundConnectivityParams.poolingOptimizationParamsStruct.poolingModel, ...
    'surroundConnectableConeTypes', rfSurroundConnectivityParams.poolingOptimizationParamsStruct.coneTypesToBePooled, ...
    'cachedData', RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.cache(...
                    theMRGCMosaic, theTargetRGCindex, ...
                    rfSurroundConnectivityParams.poolingOptimizationParamsStruct) ...
    );

  % Assemble modelVariables struct
  modelVariables = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.modelVariables(...
    modelConstants, anatomicalRcDegs, ...
    rfSurroundConnectivityParams.poolingOptimizationParamsStruct.poolingModel, ...
    userSuppliedInitialValuesForModelVariables);

  rmseSequence = [];
  poolingParamsValueSequence = [];

  % Axes handles to display the fitted model params
  hFig = figure(555); clf;
  set(hFig, 'Position', [750 700 2000 900], 'Color', [1 1 1]);

  axDoGFitToCompositeSTFparams = subplot(2,4,1);
  axSurroundWeightsParams = subplot(2,4,2);
  axSTF = subplot(2,4,3);
  axWeightsMap = subplot(2,4,4);
  axPoolingParamsSequence = subplot(2,4,5);
  axRMSEsequence = subplot(2,4,[6 7 8]);

  % Compute initial visual STF to make sure everything is running ok
  [theInitialRMSE, theInitialSTFdata, theInitialPooledConeIndicesAndWeights] = theObjectiveFunction(modelVariables.initialValues);

  % Visualize initial setup
  RGCMosaicConstructor.visualize.currentOptimizationStage(...
    theMRGCMosaic, theInitialPooledConeIndicesAndWeights, theInitialSTFdata, ...
    modelConstants.cachedData.RFcenterPosDegs, modelConstants.cachedData.maxSurroundSupportDegs, ...
    axSTF, axWeightsMap, 'initial');


  % Create optimization problem
  problem = createOptimProblem('fmincon',...
    'objective', @theObjectiveFunction, ...
    'x0', modelVariables.initialValues, ...
    'lb', modelVariables.lowerBounds, ...
    'ub', modelVariables.upperBounds, ...
    'options', optimoptions(...
      'fmincon',...
      'Display', rfSurroundConnectivityParams.fitParams.weightsOptimization.displayProgress,...
      'Algorithm','interior-point', ... %'interior-point', 'sqp'...
      'GradObj', 'off', ...
      'DerivativeCheck', 'off', ...
      'MaxFunEvals', 10^5, ...
      'MaxIter', rfSurroundConnectivityParams.fitParams.weightsOptimization.maxIterations, ...
      'UseParallel', rfSurroundConnectivityParams.fitParams.weightsOptimization.useParallel ...
      ) ...
    );

  % Start timer 
  tStart = cputime; 

  if (rfSurroundConnectivityParams.fitParams.weightsOptimization.multiStartsNum == 1)
    % Run the solver
    modelVariables.finalValues = fmincon(problem);
    if (isempty(modelVariables.finalValues))
          % Re-run the  solver
          fprintf(2, 'Solver failed. Trying once again\n');
          modelVariables.finalValues = fmincon(problem);
    end
  else
    ms = MultiStart(...
      'Display',  rfSurroundConnectivityParams.fitParams.weightsOptimization.displayProgress, ...
      'StartPointsToRun','bounds', ...  % run only initial points that are feasible with respect to bounds
      'UseParallel', false);
    
    % Run the multi-start solver
    modelVariables.finalValues = run(ms, problem, rfSurroundConnectivityParams.fitParams.weightsOptimization.multiStartsNum); 

    if (isempty(modelVariables.finalValues))
        % Re-run the multi-start solver
        fprintf(2, 'Solver failed. Trying once again\n');
        modelVariables.finalValues = run(ms, problem, rfSurroundConnectivityParams.fitParams.weightsOptimization.multiStartsNum);
    end
  end

  % Compute final visual STF
  [theFinalRMSE, theFinalSTFdata, theFinalPooledConeIndicesAndWeights] = theObjectiveFunction(modelVariables.finalValues);

  % Visualize the final values
  RGCMosaicConstructor.visualize.fittedModelParams(axSurroundWeightsParams, modelVariables, sprintf('%s\n(final)', mfilename));

  % Report elapsed time
  message = sprintf('Finished optimizing surround pooling in %f seconds.', cputime-tStart);
  disp(message);

  % Return struct
  dataOut.modelVariables = modelVariables;
  dataOut.modelConstants = modelConstants;
  dataOut.theOptimalOrientation = theOptimalOrientation;
  dataOut.visualRcDegs = visualRcDegs;
  dataOut.anatomicalRcDegs = anatomicalRcDegs;
  dataOut.theFinalRMSE = theFinalRMSE;
  dataOut.theFinalSTFdata = theFinalSTFdata;
  dataOut.theFinalPooledConeIndicesAndWeights = theFinalPooledConeIndicesAndWeights;



  %  ------- theObjectiveFunction --------
  function [theCurrentRMSE, theCurrentSTFdata, theCurrentPooledConeIndicesAndWeights] = ...
            theObjectiveFunction(theCurrentPoolingParamValues)

      % Compute current pooled cone indices and their weights
      theCurrentPooledConeIndicesAndWeights = modelConstants.poolingModel.weightsComputeHandle(modelConstants, theCurrentPoolingParamValues);

      fixedOptimalOrientation = theOptimalOrientation;   % fix it to the value determined by the center STF
      % fixedOptimalOrientation = [];  % determine it all the time from the current composite STF

      useEntireSTF = true;
      if (useEntireSTF)
        deltaThresholdToDeclareLocalMinInSTF = inf;
      else
        deltaThresholdToDeclareLocalMinInSTF = rfSurroundConnectivityParams.fitParams.DifferenceOfGaussiansToCompositeSTF.deltaThresholdToDeclareLocalMinInSTF;
      end

      
      theCurrentSTFdata = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.currentSTF(...
          theCurrentPooledConeIndicesAndWeights, inputConeMosaicSTFresponses, visualRcDegs, fixedOptimalOrientation, ...
          targetVisualSTFparams, ...
          rfSurroundConnectivityParams.fitParams.DifferenceOfGaussiansToCompositeSTF.targetVisualSTFparamFractionalTolerances,  ...
          rfSurroundConnectivityParams.fitParams.DifferenceOfGaussiansToCompositeSTF.multiStartsNum, ...
          deltaThresholdToDeclareLocalMinInSTF, ...
          stimParams, ...
          axDoGFitToCompositeSTFparams);

      % Compute RMSE
      RsRcRatioResidual = ...
        theCurrentSTFdata.fittedDoGModelRsRcRatio / targetVisualSTFparams.surroundToCenterRcRatio - 1;
      SCintSensRatioResidual = ...
        theCurrentSTFdata.fittedDoGModelSCIntSensRatio / targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio - 1;

      residualWeights = [0 0  1];
      allResiduals = [RsRcRatioResidual SCintSensRatioResidual theCurrentSTFdata.residual];
      theCurrentRMSE = norm(allResiduals.*residualWeights);
      rmseSequence(size(rmseSequence,1)+1,:) = cat(2, theCurrentRMSE, allResiduals);

      if (rmseSequence(end,3) == min(rmseSequence(:,3)))
        poolingParamsValueSequence(size(poolingParamsValueSequence,1)+1,:) = theCurrentPoolingParamValues([3 5]);
        scatter(axPoolingParamsSequence, poolingParamsValueSequence(:,1), poolingParamsValueSequence(:,2), 100, ...
          'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4);
        xlabel(axPoolingParamsSequence, 'V ratio');
        ylabel(axPoolingParamsSequence, 'R ratio');
        title('H1 model params')
        set(axPoolingParamsSequence, 'XLim', [0 1.2], 'YLim', [0 0.4], 'FontSize', 12);
      end


      % Update plots
      modelVariables.finalValues = theCurrentPoolingParamValues;
      if (strcmp(modelConstants.poolingModel.name, 'PackerDacey2002H1FixedCellIndex'))
        plotTitle = sprintf('%s (H1-%d)', modelConstants.poolingModel.name, modelConstants.poolingModel.fixedH1CellIndex);
      else
        plotTitle = sprintf('%s', modelConstants.poolingModel.name);
      end

      RGCMosaicConstructor.visualize.fittedModelParams(axSurroundWeightsParams, modelVariables, plotTitle );


      currentIteration = size(rmseSequence,1);
      plot(axRMSEsequence, 1:currentIteration, 0*rmseSequence(:,1), 'k--', 'LineWidth', 1.0);
      hold(axRMSEsequence, 'on');
      p1 = plot(axRMSEsequence, 1:currentIteration, rmseSequence(:,1), 'k-', 'LineWidth', 2.0);
      p2 = plot(axRMSEsequence, 1:currentIteration, rmseSequence(:,2), 'ro',  'MarkerSize', 4, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0);
      p3 = plot(axRMSEsequence, 1:currentIteration, rmseSequence(:,3), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', [0 1 1], 'LineWidth', 1.0);
      p4 = plot(axRMSEsequence, 1:currentIteration, rmseSequence(:,4), 'k:', 'MarkerSize', 4, 'MarkerFaceColor', [0 1 1], 'LineWidth', 1.0);
      hold(axRMSEsequence, 'off');
      legend(axRMSEsequence, [p1 p2 p3 p4], {'total residual', 'RsRc residual', 'SCintSens residual', 'DoGfit-STF residual'});
      xlabel(axRMSEsequence, 'iteration no');
      ylabel(axRMSEsequence, 'residuals');
      title(axRMSEsequence,sprintf('Current residual: %f', rmseSequence(end,1)));
      %iterationsToConsider = max([1 currentIteration-20]):currentIteration;
      %yLimRange = max([1.0 min([10 2*mean(mean(abs(rmseSequence(iterationsToConsider,:))))])]);
      yLimRange = 0.3;

      if (yLimRange > 10)
        deltaTick = 2;
      elseif (yLimRange > 5)
        deltaTick = 1.0;
      elseif (yLimRange > 2)
        deltaTick = 0.5;
      elseif (yLimRange > 1)
        deltaTick = 0.2;
      elseif (yLimRange > 0.5)
        deltaTick = 0.1;
      else
        deltaTick = 0.05;
      end

      set(axRMSEsequence, 'YTick', -100:deltaTick:100, 'YLim', yLimRange*[-1 1], ...
        'XLim', [1 max(10,currentIteration)], 'FontSize', 13);
      grid(axRMSEsequence, 'on');

      if (theCurrentRMSE == min(squeeze(rmseSequence(:,1))))
          RGCMosaicConstructor.visualize.currentOptimizationStage(...
            theMRGCMosaic, theCurrentPooledConeIndicesAndWeights, theCurrentSTFdata, ...
            modelConstants.cachedData.RFcenterPosDegs, modelConstants.cachedData.maxSurroundSupportDegs, ...
            axSTF, axWeightsMap, sprintf('iter: %d (total residual:%2.3f)', currentIteration,theCurrentRMSE));
      end

      drawnow;

  end %  ------- theObjectiveFunction --------

end