function [RcDegs, visualRFcenterConeMap, retinalRFcenterConeMap, anatomicalConeCharacteristicRadiusDegs] = estimateVisualRcFromNumberOfConesInRFcenter(cm, conesNumPooledByTheRFcenter, theCircularPSFData)

        % Compute the visual Rc based on the number of cones in the RF center 
        % their spacing and the PSF, all for the current eccentricity

        % Sort cones according to their distance from the mosaic center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');

        % Estimate mean anatomical cone aperture in the mosaic'c center
        sourceConesIndices = idx(1:6);
        meanConeApertureDegsInMosaicCenter = mean(cm.coneApertureDiametersDegs(sourceConesIndices));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * meanConeApertureDegsInMosaicCenter;

        % Compute the retinal RF center cone map
        rfCenterPooledConeIndices = idx(1:conesNumPooledByTheRFcenter);
        meanRFCenterConePos = mean(cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:),1);
        conePosDegsRelativeToCenter = bsxfun(@minus, cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:), meanRFCenterConePos);
        
        spatialSupportDegs = [theCircularPSFData.supportX(:) theCircularPSFData.supportY(:)]/60;
        [~, retinalRFcenterConeMap] = RetinaToVisualFieldTransformer.computeRetinalRFRcDegsFromItsPooledConeInputs(...
            anatomicalConeCharacteristicRadiusDegs, conePosDegsRelativeToCenter, spatialSupportDegs);

        % Convolve the retinal RF center cone map with the PSF
        visualRFcenterConeMap = conv2(theCircularPSFData.data, retinalRFcenterConeMap, 'same');

        % Fit a 2D Gaussian ellipsoid to the visually projected RF center cone map and extract
        % the characteristic radii of that Gaussian ellipsoid
        [~,visualConeCharacteristicMinorMajorRadiiDegs] = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            theCircularPSFData.supportX, theCircularPSFData.supportY, visualRFcenterConeMap);

        % Use the mean of the radii of the Gaussian ellipsoid as the visual RF Rc
        RcDegs = mean(visualConeCharacteristicMinorMajorRadiiDegs);
end
