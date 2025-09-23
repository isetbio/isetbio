function currentOptimizationStage(theMRGCMosaic, currentPooledConeIndicesAndWeights, currentSTFdata, ...
  RFcenterPosDegs, maxSurroundSupportDegs, axSTF, axWeightsMap, plotTitle)

  plot(axSTF,currentSTFdata.spatialFrequencySupport,  currentSTFdata.visualSTF, 'ko', 'MarkerSize', 16, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 1.0);
  hold(axSTF, 'on');
  plot(axSTF, currentSTFdata.spatialFrequencySupportHR,  currentSTFdata.visualCompositeSTFHR, 'k-', 'LineWidth', 1.5);
  plot(axSTF, currentSTFdata.spatialFrequencySupportHR,  currentSTFdata.visualCenterSTFHR, 'r-','LineWidth', 1.5);
  plot(axSTF, currentSTFdata.spatialFrequencySupportHR,  currentSTFdata.visualSurroundSTFHR, 'b-', 'LineWidth', 1.5);
  set(axSTF,  'XScale', 'log');
  axis(axSTF, 'square')
  hold(axSTF, 'off');
  title(axSTF, sprintf('Rs/Rc = %2.2f, S/C intS = %2.2f',...
       currentSTFdata.fittedDoGModelRsRcRatio, currentSTFdata.fittedDoGModelSCIntSensRatio ));
  
  currentSurroundWeights = currentPooledConeIndicesAndWeights.surroundConeWeights;
  idxSurround = find(currentSurroundWeights > max(currentSurroundWeights(:))*3/100);
  surroundConesNum = numel(idxSurround);

  centerWeights = currentPooledConeIndicesAndWeights.centerConeWeights;
  idxCenter = find(centerWeights >= mRGCMosaic.sensitivityAtPointOfOverlap);
  centerConesNum = numel(idxCenter);


  plot(axWeightsMap, ...
    theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(currentPooledConeIndicesAndWeights.surroundConeIndices(idxSurround),1), ...
    theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(currentPooledConeIndicesAndWeights.surroundConeIndices(idxSurround),2), ...
    'bo', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.8 1]);
  hold(axWeightsMap, 'on');
  plot(axWeightsMap, ...
    theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(currentPooledConeIndicesAndWeights.centerConeIndices(idxCenter),1), ...
    theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(currentPooledConeIndicesAndWeights.centerConeIndices(idxCenter),2), ...
    'ro', 'LineWidth', 1.0, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
  axis(axWeightsMap, 'square');

  set(axWeightsMap, ...
    'XLim', RFcenterPosDegs(1) + 0.9*maxSurroundSupportDegs*[-1 1], ...
    'YLim', RFcenterPosDegs(2) + 0.9*maxSurroundSupportDegs*[-1 1]);

  plotTitle = sprintf('%s (center cones: %d, surround cones: %d)', plotTitle, centerConesNum, surroundConesNum);
  hold(axWeightsMap, 'off');
  title(axWeightsMap, plotTitle);
end
