function synthesizedRFParams = synthesizeRFparams(conePositionsMicrons, coneSpacingsMicrons, midgetRGCconnectionMatrix, ...
    deconvolutionOpticsParams)
    
    % Determine rf center size and rf center position from the connectivity matrix
    [rfCenterRadiiMicrons, rfCenterPositionsMicrons] = ...
        computeRFcenterSizeAndPositionsFromConnectivityMatrix(...
            conePositionsMicrons, coneSpacingsMicrons, ...
            midgetRGCconnectionMatrix);
        
    rfCenterRadiiMicrons = mean(rfCenterRadiiMicrons,2);
    %rfCenterRadiiMicrons = min(rfCenterRadiiMicrons,[],2);
    rfCenterRadiiMicrons = max(rfCenterRadiiMicrons,[],2);
    
    % Compute visual and retinal RF params from the retinal rf center radii
    % and positions (both in microns)
    ck = CronerKaplanRGCModel('generateAllFigures', false, 'instantiatePlotLab', false);
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(...
        rfCenterRadiiMicrons, rfCenterPositionsMicrons, deconvolutionOpticsParams);
    
    % position as computed by the summed inputs to the center
    synthesizedRFParams.centerPositionMicrons = rfCenterPositionsMicrons;  
end

function [rfCenterRadiiMicrons, rfCenterPositionsMicrons] = ...
    computeRFcenterSizeAndPositionsFromConnectivityMatrix(...
            conePositionsMicrons, coneSpacingsMicrons, ...
            midgetRGCconnectionMatrix)
        
    % Preallocate memory
    rgcsNum = size(midgetRGCconnectionMatrix,2);
    rfCenterRadiiMicrons = zeros(rgcsNum,2);
    rfCenterPositionsMicrons = zeros(rgcsNum,2);
    
    deltaX = min(coneSpacingsMicrons)/20;
    
    % Compute RGC RF center positions and radii
    for RGCindex = 1:rgcsNum
        
        connectivityVector = full(squeeze(midgetRGCconnectionMatrix(:, RGCindex)));
        inputIDs = find(connectivityVector>0);
        inputsNum = numel(inputIDs);
        if (inputsNum == 0)
            error('RGC %d has zero inputs\n', RGCindex);
        end
        
        % Compute RF center position and size based on the connectivity matrix
        switch inputsNum
            case 1
                 rfCenterPosition = conePositionsMicrons(inputIDs,:);
                 rfSigmas = 0.5*coneSpacingsMicrons(inputIDs)/3*[1 1];
                 
            case 2
                [rfCenterPosition, rfSigmas] = computeRFcenterAndSizeFor2inputRGCFromConnectivityMatrix(...
                    connectivityVector(inputIDs), conePositionsMicrons(inputIDs,:), coneSpacingsMicrons(inputIDs));
                
            otherwise
                [rfCenterPosition, rfSigmas] = computeRFcenterAndSizeFromConnectivityMatrixViaGaussianFitting(...
                    connectivityVector(inputIDs), conePositionsMicrons(inputIDs,:), coneSpacingsMicrons(inputIDs), deltaX);
        end
        
        fprintf('RF center size for RGC %d (%d-input):%2.2f um\n', RGCindex, inputsNum, mean(rfSigmas));
        
        rfCenterPositionsMicrons(RGCindex,:) = rfCenterPosition;
        rfCenterRadiiMicrons(RGCindex,:) = rfSigmas;
    end 
end

function [rfCenterPosition, rfSigmas] = computeRFcenterAndSizeFor2inputRGCFromConnectivityMatrix(...
            connectivityVector, conePositionsMicrons, coneSpacingsMicrons)
        
       if (numel(connectivityVector) ~= 2)
           error('Must have 2 only inputs');
       end
       
       rfCenterPosition = [];
       for coneIndex = 1:numel(connectivityVector)
           if (isempty(rfCenterPosition))
               rfCenterPosition = conePositionsMicrons(coneIndex,:) * connectivityVector(coneIndex);
           else
               rfCenterPosition = rfCenterPosition + conePositionsMicrons(coneIndex,:) * connectivityVector(coneIndex);
           end
           
       end
       coneApertureRadius = 0.7*0.5*coneSpacingsMicrons(1);
       coneSigma = coneApertureRadius/3;
       rfCenterPosition = sum(rfCenterPosition,1)/sum(connectivityVector);
       rfSigmas = [coneSigma coneSpacingsMicrons(1)*0.5+coneSigma];

end


function [rfCenterPosition, rfSigmas] = computeRFcenterAndSizeFromConnectivityMatrixViaGaussianFitting(connectivityVector, conePositionsMicrons, coneSpacingsMicrons, deltaX)
    % Determine spatial support
    minXY = min(conePositionsMicrons,[],1);
    maxXY = max(conePositionsMicrons,[],1);
    s = max(coneSpacingsMicrons);
    xMin = minXY(1)-s; xMax = maxXY(1)+s;
    yMin = minXY(2)-s; yMax = maxXY(2)+s;
    xAxis = (xMin: deltaX: xMax);
    yAxis = (yMin: deltaX: yMax);
    [X,Y] = meshgrid(xAxis,yAxis);
    theRF = [];
    
    for coneIndex = 1:numel(connectivityVector)
        cP = squeeze(conePositionsMicrons(coneIndex,:));
        coneApertureRadius = 0.7*0.5*coneSpacingsMicrons(coneIndex);
        coneSigma = coneApertureRadius/3;
        coneProfile = (exp(-((X-cP(1))/coneSigma).^2) .* exp(-((Y-cP(2))/coneSigma).^2));
        
        % Flat-top sensitivity
        %coneProfile(coneProfile>=exp(-1)*exp(-1)) = exp(-1)*exp(-1);
        %coneProfile = coneProfile / max(coneProfile(:));
        
        if (isempty(theRF))
            theRF = coneProfile * connectivityVector(coneIndex);
        else
            theRF = theRF + coneProfile * connectivityVector(coneIndex);
        end 
    end
    
    % Fit a 2D Gaussian with different minor/major axes
    center = mean(conePositionsMicrons,1);
    [rfCenterPosition, rfSigmas] = fitElliptical2DGausianToRF(X, Y, theRF/max(theRF(:)), deltaX, center);
    
end


