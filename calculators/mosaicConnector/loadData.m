function [RGCRFPositions, RGCRFSpacings, conePositions, conesToRGCratios] = ...
        loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, eccentricitySamplesNumRGC, ...
        maxMovementPercentileCones, maxMovementPercentileRGC,  bestIterationCones,  bestIterationRGC, analyzedRadiusDeg)
    
    micronsPerDeg = 300;
    workingRadiusMicrons = analyzedRadiusDeg*micronsPerDeg;
    
    % Load the positions of the mRGC mosaic
    neuronalType = 'mRGC';
    RGCRFPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumRGC, maxMovementPercentileRGC,  bestIterationRGC);
    
    % Load the positions of the cone mosaic
    neuronalType = 'cone';
    conePositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, maxMovementPercentileCones,  bestIterationCones);
  
    % Compute spacing and cone-to-RGC ratios for each  RGCRF position
    [RGCRFSpacings, conesToRGCratios] = mRGCStats(RGCRFPositions, 128, whichEye);
      
    % Only keep cones  within the working radius
    eccCones = sqrt(sum(conePositions.^2,2));
    idx = find((eccCones<workingRadiusMicrons) & (conePositions(:,1)>=-20));
    conePositions = conePositions(idx,:);
    
    % Only keep RGC within the working radius
    eccRGC = sqrt(sum(RGCRFPositions.^2,2));
    idx = find((eccRGC<workingRadiusMicrons) & (RGCRFPositions(:,1)>=-20));
    RGCRFPositions = RGCRFPositions(idx,:);
    RGCRFSpacings = RGCRFSpacings(idx);
    conesToRGCratios = conesToRGCratios(idx);
    
    minRGCEccDegs = min(RGCRFPositions(:,1))/micronsPerDeg;
    maxRGCEccDegs = max(RGCRFPositions(:,1))/micronsPerDeg;
    
    minRGCEccCones = min(conePositions(:,1))/micronsPerDeg;
    maxRGCEccCones = max(conePositions(:,1))/micronsPerDeg;
    
    fprintf('Loaded %2.0f RGCs from %2.2f-%2.2f degs\n', size(RGCRFPositions,1), minRGCEccDegs, maxRGCEccDegs);
    fprintf('Loaded %2.0f cones from %2.2f-%2.2f degs\n', size(conePositions,1), minRGCEccCones, maxRGCEccCones);    
end

function rfPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum, maxMovementPercentile, bestIteration)
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    
    mosaicFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_maxMovPrctile%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, maxMovementPercentile));

    load(mosaicFileName, 'rfPositionsHistory', 'reTriangulationIterations');
    if (~isinf(bestIteration))
        [~, targetIterationIndex] = min(abs(reTriangulationIterations - bestIteration));
    else
        targetIterationIndex = numel(reTriangulationIterations);
    end
    
    fprintf('Loading %s mosaic with %d neurons from iteration %d\n', neuronalType, size(rfPositionsHistory,2), reTriangulationIterations (targetIterationIndex));
    rfPositions = double(squeeze(rfPositionsHistory(targetIterationIndex,:,:)));   
end