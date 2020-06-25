% Phase 6: Compute weights of cone inputs to mRGC  RF centers and surrounds
function runPhase6(runParams)

    % Connectivity datafile
    connectivityFile = fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile));
    
    % Compute RF center radius and position from the connectivity matrix (both in microns) for each mRGC
    % in the passed patch
    patchEccDegs = runParams.patchEccDegs;
    patchSizeDegs = runParams.patchSizeDegs;
    [rfCenterRadiiMicrons, rfCenterPositionsMicrons, rgcIndices] = ...
        computeRFcenterSizeAndPositionsFromConnectivityMatrix(connectivityFile, patchEccDegs, patchSizeDegs);
    
    % Compute visual and retinal RF params from the retinal rf center radii
    % and positions (both in microns)
    ck = CronerKaplanRGCModel('generateAllFigures', false, 'instantiatePlotLab', false);
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(rfCenterRadiiMicrons, rfCenterPositionsMicrons);
    
    % Add the corresponding rgcIndices, and patch ecc/size
    synthesizedRFParams.rgcIndices = rgcIndices;     % indices of RGCs within the target patch
    synthesizedRFParams.centerPositionMicrons = rfCenterPositionsMicrons;  % position as computed by the summed inputs to the center
    synthesizedRFParams.patchEccDegs = patchEccDegs;
    synthesizedRFParams.patchSizeDegs = patchSizeDegs;
       
    % Compute weights of cone inputs to the RF center and to the RF
    % surround based on the midgetRGCconnectionMatrix and the
    % synthesizedRFparams
    load(connectivityFile , ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios', 'midgetRGCconnectionMatrix');
        
    [midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround] = ...
        computeWeightedConeInputsToRGCCenterSurroundSubregions(...
            synthesizedRFParams.retinal, ...
            synthesizedRFParams.centerPositionMicrons, ...
            synthesizedRFParams.eccDegs, ...
            synthesizedRFParams.rgcIndices, ...
            conePositionsMicrons, coneTypes, midgetRGCconnectionMatrix);
    
    % Save synthesizedRF params (cone weights to center/surround regions)
    saveDataFile = sprintf('%s_inPatchAt_%2.1f_%2.1fdegs_WithSize_%2.2f_%2.2f.mat', ...
        runParams.outputFile, patchEccDegs(1), patchEccDegs(2), patchSizeDegs(1), patchSizeDegs(2));
        
    save(fullfile(runParams.outputDir, saveDataFile), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', 'desiredConesToRGCratios', ...
            'midgetRGCconnectionMatrixCenter', 'midgetRGCconnectionMatrixSurround', ...
            'synthesizedRFParams', 'patchEccDegs', 'patchSizeDegs', '-v7.3');
end

