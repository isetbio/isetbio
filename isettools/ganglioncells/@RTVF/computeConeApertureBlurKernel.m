function computeConeApertureBlurKernel(obj)

    theConeIndices = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;

    % Find the blurZone in which theConeIndices belong to
    targetZoneIndex = coneApertureBlurZone(obj.coneMosaic, theConeIndices);

    % Retrieve the blurApertureDiameter for the target zone
    blurApertureDiameterDegs = obj.coneMosaic.blurApertureDiameterDegsZones(targetZoneIndex);

    % Compute the aperture blur kernel
    lowOpticalImageResolutionWarning = true;
    oiResDegs = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs(2) - ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs(1);
    obj.coneApertureBlurKernel = obj.coneMosaic.generateApertureKernel(...
        blurApertureDiameterDegs(1), oiResDegs, lowOpticalImageResolutionWarning);

end

function targetZoneIndex = coneApertureBlurZone(theConeMosaic, theConeIndices)
    % Find the blur aperture diameter zone for theConeIndices 
    coneApertureBlurZonesNum = numel(theConeMosaic.blurApertureDiameterMicronsZones);
    inputConesNumInZone = zeros(1,coneApertureBlurZonesNum);

    for zoneIndex = 1:coneApertureBlurZonesNum 
        % Determine which cones should receive this blur.
        coneIDsInZone = theConeMosaic.coneIndicesInZones{zoneIndex};
        inputConesNumInZone(zoneIndex) = numel(find(ismember(theConeIndices,coneIDsInZone) == 1));
    end

    [~, targetZoneIndex] = max(inputConesNumInZone);
end
