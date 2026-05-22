function mergedConeMosaicActivation = computeForOIensemble(obj, oiEnsemble, theMergingWeights, theScene)

    % Currently only works for static images without fixational EMs
    opticsSamplingPositionsNum = numel(oiEnsemble);

    for oiPos = 1:opticsSamplingPositionsNum
        % Retrieve the OI at the current position
        theOI = oiEnsemble{oiPos};

        % Compute the retinal image of the entire scene under this OI
        theRetinalImage = oiCompute(theOI, theScene, 'pad value','mean');

        % Compute the cone mosaic activation for the current OI
        theConeMosaicActivation = obj.compute(...
            theRetinalImage, ...
            'opticalImagePositionDegs', [0 0]);

        % Merge with the cone mosaic activations for the other OI positions
        if (oiPos == 1)
            mergedConeMosaicActivation = theMergingWeights(oiPos,:) .* (theConeMosaicActivation(:))';
        else
            mergedConeMosaicActivation =  mergedConeMosaicActivation + ...
                theMergingWeights(oiPos,:) .* (theConeMosaicActivation(:))';
        end

    end % oiPos

    mergedConeMosaicActivation = reshape(mergedConeMosaicActivation, size(theConeMosaicActivation));
end