function [RcDegs, rfRotationDegs, flatTopGaussianExponent, ...
    visualRFcenterConeMap, retinalRFcenterConeMap, ...
    anatomicalConeCharacteristicRadiusDegs] = estimateVisualRcFromNumberOfConesInRFcenter(cm, conesNumPooledByTheRFcenter, theCircularPSFData, spatialSupportDegs)

        % Compute the visual Rc based on the number of cones in the RF center 
        % their spacing and the PSF, all for the current eccentricity

        % Sort cones according to their distance from the mosaic center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');

        % Estimate mean anatomical cone aperture in the mosaic'c center
        sourceConesIndices = idx(1:6);
        maxConeApertureDegsInMosaicCenter = max(cm.coneApertureDiametersDegs(sourceConesIndices));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * maxConeApertureDegsInMosaicCenter;

        % Compute the retinal RF center cone map
        rfCenterPooledConeIndices = idx(1:conesNumPooledByTheRFcenter);
        meanRFCenterConePos = mean(cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:),1);
        conePosDegsRelativeToCenter = bsxfun(@minus, cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:), meanRFCenterConePos);    
        [~, retinalRFcenterConeMap] = RetinaToVisualFieldTransformer.computeRetinalRFRcDegsFromItsPooledConeInputs(...
                anatomicalConeCharacteristicRadiusDegs, conePosDegsRelativeToCenter, spatialSupportDegs);


        % Convolve the retinal RF center cone map with the PSF
        visualRFcenterConeMap = conv2(retinalRFcenterConeMap, theCircularPSFData.data, 'same');

        % Fit a 2D Gaussian ellipsoid to the visually projected RF center cone map and extract
        % the characteristic radii of that Gaussian ellipsoid
        rfSupportX = spatialSupportDegs(:,1);
        rfSupportY = spatialSupportDegs(:,2);
        [~,visualConeCharacteristicMinorMajorRadiiDegs, rfRotationDegs, flatTopGaussianExponent] = ...
            RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                rfSupportX, rfSupportY, visualRFcenterConeMap, ...
                'flatTopGaussian', true);

        
        RcDegs = visualConeCharacteristicMinorMajorRadiiDegs;
end
