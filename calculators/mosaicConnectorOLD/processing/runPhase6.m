% Phase 6: Compute weights of cone inputs to mRGC  RF centers and surrounds
function runPhase6(runParams)

    connectivityFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
    load(connectivityFile, ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios', 'midgetRGCconnectionMatrix');
    
    
    deconvolutionOpticsParams = runParams.deconvolutionOpticsParams;
    coneMosaicEccMicrons = 1e3*WatsonRGCModel.rhoDegsToMMs(runParams.patchEccDegs);
    coneMosaicSizeMicrons = WatsonRGCModel.sizeDegsToSizeRetinalMicrons(runParams.patchSizeDegs, runParams.patchEccDegs);
    
    [midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, ...
        synthesizedRFParams] = ...
        computeWeightedConeInputsToRGCCenterSurroundSubregions(...
            conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            midgetRGCconnectionMatrix, ...
            coneMosaicEccMicrons, coneMosaicSizeMicrons, ...
            deconvolutionOpticsParams);
    
   
    % Assemble cone weights data file name
    qNum = numel(runParams.deconvolutionOpticsParams.quadrantsToAverage);
    quadrantsToAverage = '';
    for qIndex = 1:qNum
        quadrantsToAverage = sprintf('_%s', quadrantsToAverage, runParams.deconvolutionOpticsParams.quadrantsToAverage{qIndex});
    end

    saveDataFile = sprintf('%s_inPatchAt_%2.1f_%2.1fdegs_WithSize_%2.2f_%2.2f_ForSubject_%d_AndQuadrants%s.mat', ...
        runParams.outputFile, runParams.patchEccDegs(1), runParams.patchEccDegs(2), ...
        runParams.patchSizeDegs(1), runParams.patchSizeDegs(2), ...
        runParams.deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage, ...       // Deconvolution model: which subject
        quadrantsToAverage...                                            // Deconvolution model: which quadrant to use/average
      );
        
    patchEccDegs = runParams.patchSizeDegs;
    patchSizeDegs = runParams.patchEccDegs;
    PolansSubjectIDsAveraged = runParams.deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage;
    quadrantsAveraged = runParams.deconvolutionOpticsParams.quadrantsToAverage;
    
     % Save cone weights to center/surround regions
    save(fullfile(runParams.outputDir, saveDataFile), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios', ...
            'midgetRGCconnectionMatrixCenter', 'midgetRGCconnectionMatrixSurround', ...
            'synthesizedRFParams', 'patchEccDegs', 'patchSizeDegs', ...
            'PolansSubjectIDsAveraged', 'quadrantsAveraged', '-v7.3');
end
