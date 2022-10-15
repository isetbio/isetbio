function retinalRFcenterMaps = computeRetinalRFcenterMaps(obj, marginDegs, spatialSupportSamplesNum)

    % Preallocate memory
    mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    retinalRFcenterMaps = cell(1, mRGCsNum);

    parfor iRGC = 1:mRGCsNum
        % Find this RGC's center input cone indices and their weights
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        inputConeIndices = find(connectivityVector > 0.0001);

        %if (all(connectivityVector(inputConeIndices)==1))
            % Add some nearby cone positions with zero weights to help with
            % fitting purposes
            centroid = mean(obj.inputConeMosaic.coneRFpositionsDegs(inputConeIndices,:),1);
            [~, otherConeIndices] = MosaicConnector.pdist2(obj.inputConeMosaic.coneRFpositionsDegs, centroid, ...
                'smallest', numel(inputConeIndices)+max([7 ceil(2*pi*sqrt(numel(inputConeIndices)))]) );
            additionalConeIndices = setdiff(otherConeIndices, inputConeIndices);
            inputConeIndices = [inputConeIndices; additionalConeIndices];

       % end
        allInputConeWeights = connectivityVector(inputConeIndices);

        % Center input cone positions and characteristic radii
        allInputConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(inputConeIndices,:);
        allInputConeRcDegs = obj.inputConeMosaic.coneApertureDiametersDegs(inputConeIndices) * ...
                             obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;

        % Determine x and y range based on this mRGC center cone inputs
        [spatialSupportDegsX, spatialSupportDegsY] = spatialSupportForCurrentRGC(allInputConePositionsDegs, marginDegs, spatialSupportSamplesNum);
        [X,Y] = meshgrid(spatialSupportDegsX, spatialSupportDegsY);
        
        % Assemble the center RF map from its cone RF maps
        theRFcenterMap = X*0;
        coneApertures = zeros(numel(inputConeIndices), size(X,1), size(X,2));
        for iInput = 1:numel(inputConeIndices)
             inputConePosDegs = allInputConePositionsDegs(iInput,:);
             inputConeRcDegs = allInputConeRcDegs(iInput);
             XX = X - inputConePosDegs(1);
             YY = Y - inputConePosDegs(2);
             R = sqrt(XX.^2 + YY.^2);
             coneApertures(iInput,:,:) = allInputConeWeights(iInput)*exp(-(R/inputConeRcDegs).^2);
             theRFcenterMap = theRFcenterMap + squeeze(coneApertures(iInput,:,:));
        end 

        % Save the RF map
        retinalRFcenterMaps{iRGC} = struct(...
             'centerRF', theRFcenterMap, ...
             'coneApertures', coneApertures, ...
             'inputConeIndices', inputConeIndices, ...
             'inputConeWeights', allInputConeWeights, ...
             'spatialSupportDegsX', spatialSupportDegsX, ...
             'spatialSupportDegsY', spatialSupportDegsY...
             );
    end % iRGC
end


function [spatialSupportDegsX, spatialSupportDegsY] = spatialSupportForCurrentRGC(allInputConePositionsDegs, marginDegs, spatialSamplesNum)
    minXY = min(allInputConePositionsDegs,[],1);
    maxXY = max(allInputConePositionsDegs,[],1);
    meanXY = mean(allInputConePositionsDegs,1);
    xRange = maxXY(1)-minXY(1);
    yRange = maxXY(2)-minXY(2);
    xyRange = max([xRange yRange]);
    x1 = meanXY(1)-xyRange/2-1.2*marginDegs;
    x2 = meanXY(1)+xyRange/2+1.2*marginDegs;
    y1 = meanXY(2)-xyRange/2-1.2*marginDegs;
    y2 = meanXY(2)+xyRange/2+1.2*marginDegs;

    % Compute spatial support
    spatialSupportDegsX = linspace(x1,x2, spatialSamplesNum);
    spatialSupportDegsY = linspace(y1,y2, spatialSamplesNum);
end