function coneDensityMMs2 = correctionForISETBioFovealConeDensity(coneDensityMMs2, eccDegs)
% Apply correction for the fact that the isetbio max cone density (18,800 cones/deg^2) 
% does not agree with Watson's (obj.dc0 =  14,804.6 cones/deg^2), and the fact that if we do not
% apply this correction we get less than 2 mRGCs/cone at foveal eccentricities. We
% apply this correction only for ecc <= 0.18 degs

    eccCorrectionLimit = 0.18;
    if (min(eccDegs) > eccCorrectionLimit)
        return;
    end
    
    % Find max cone density in ISETBio
    [~,~,ISETBioMaxConeDensityPerMMs2] = coneSizeReadData('eccentricity', 0, 'angle', 0);
    
    % Convert to density per deg^2
    [ISETBioMaxConeDensityPerDegs2, alpha] = RGCmodels.Watson.convert.densityMMs2ToDegs2(ISETBioMaxConeDensityPerMMs2, 0.0);
    
    % Determine correction factor
    correctionFactorMax = ISETBioMaxConeDensityPerDegs2 - RGCmodels.Watson.constants.dc0;
    correctionFactorMax = correctionFactorMax / alpha;
        
    % Apply corrections
    idx = find(abs(eccDegs) <= eccCorrectionLimit);
    if (~isempty(idx))
        indicesToBeCorrected = idx;
        correctionFactors = correctionFactorMax .* (eccCorrectionLimit-eccDegs(indicesToBeCorrected))/eccCorrectionLimit;
        coneDensityMMs2(indicesToBeCorrected) = coneDensityMMs2(indicesToBeCorrected) - correctionFactors;
    end
end

