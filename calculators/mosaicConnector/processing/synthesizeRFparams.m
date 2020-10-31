function synthesizedRFParams = synthesizeRFparams(conePositionsMicrons, midgetRGCconnectionMatrix, deconvolutionOpticsParams)
   
    rgcsNum = size(midgetRGCconnectionMatrix,2);
    rfCenterPositionsMicrons = zeros(rgcsNum,2);
    rfCenterInputConesNum = zeros(rgcsNum,1);
    
    parfor RGCindex = 1:rgcsNum
        % Extract connected cones
        connectivityVector = full(squeeze(midgetRGCconnectionMatrix(:, RGCindex)));
        inputIDs = find(connectivityVector>0);
        
        inputsNum = numel(inputIDs);
        if (inputsNum == 0)
            error('RGC %d has zero inputs\n', RGCindex);
        end
        
        rfCenterInputConesNum(RGCindex) = inputsNum;
        rfCenterPositionsMicrons(RGCindex,:) = ...
            sum(bsxfun(@times, conePositionsMicrons, connectivityVector),1)/sum(connectivityVector);
    end
    
    % Synthesize RF params
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(...
        rfCenterInputConesNum, rfCenterPositionsMicrons, deconvolutionOpticsParams);
    
    % position as computed by the summed inputs to the center
    synthesizedRFParams.centerPositionMicrons = rfCenterPositionsMicrons;  

end

