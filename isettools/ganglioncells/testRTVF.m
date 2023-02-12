function testRTVF
    positionDegs = [2 0]
    test1(positionDegs)
end

function test1(positionDegs)
    theConeMosaic = generateTestConeMosaic(positionDegs);

    theOpticsParams = RTVF.defaultOpticsParams;
    theOpticsParams.positionDegs = positionDegs;

    theTargetVisualRFDoGparams = RTVF.defaultTargetVisualRFDoGParams;
    
    d = sum((bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, positionDegs)).^2,2);
    [~,indexOfCenterMostCone] = min(d(:));

    theTargetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter = [indexOfCenterMostCone];
    theTargetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter = [1.0];
    
    theConeMosaic.visualize('labelConesWithIndices', theTargetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter)

    % Median value of C&K Rs/Rc ratios
    [temporalEccDegs, RcRsRatios] = RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    theTargetVisualRFDoGparams.surroundToCenterRcRatio = prctile(1./RcRsRatios, 50);
            
    % Median value of C&K S/C ratios
    [temporalEccDegs, SCratios] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
    theTargetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = prctile(SCratios,50);

    if (numel(theTargetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter) < 3)
        multiStartsNumRetinalPooling= 8;
    else
        multiStartsNumRetinalPooling= 4;
    end

    % Instantiate the default RTVF object, which consists of a cone mosaic
    % with 1 L-cone and 1 M-cone, and Polans Subject #4 optics
    myRTVF = RTVF(theConeMosaic, ...
                 theOpticsParams, ...
                 theTargetVisualRFDoGparams, ...
                 'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
                 'visualizeSpectrallyWeightedPSFs', ~true);

end


function test2()
    % Empty coneMosaic. This results in the RTVFobj generating a cone mosaic
    % with 1 L-cone (with index1) and 1 M-cone (with index 2).

    theConeMosaic = [];
    theOpticsParams = RTVF.defaultOpticsParams;
    theTargetVisualRFDoGparams = RTVF.defaultTargetVisualRFDoGParams;

    % FOR TESTING:
    % 1. Custom  optical position.
    % This affects both the PSF and the L-, M-spectral weights (which
    % differ due to changes in the macular pigment with eccentricity)
    opticalPositionDegs = [-6 0];
    theOpticsParams.positionDegs = opticalPositionDegs;

    % FOR TESTING:
    % 2. Custom  cones pooled by the RF center: 1 L-cone with unit weight
    % This affects the estimate of the visual Rc of the RF center at the
    % examined optical position
    theTargetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter = [1];
    theTargetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter = [1.0];

    % Median value of C&K Rs/Rc ratios
    [temporalEccDegs, RcRsRatios] = RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    theTargetVisualRFDoGparams.surroundToCenterRcRatio = prctile(1./RcRsRatios, 50);
            

    % Median value of C&K S/C ratios
    [temporalEccDegs, SCratios] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
    theTargetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = prctile(SCratios,50);

    % Instantiate the default RTVF object, which consists of a cone mosaic
    % with 1 L-cone and 1 M-cone, and Polans Subject #4 optics
    myRTVF = RTVF(theConeMosaic, ...
                 theOpticsParams, ...
                 theTargetVisualRFDoGparams, ...
                 'visualizeSpectrallyWeightedPSFs', ~true);

    theTargetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter = [1 2];
    theTargetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter = [1.0 1.0];

    % Instantiate the default RTVF object, which consists of a cone mosaic
    % with 1 L-cone and 1 M-cone, and Polans Subject #4 optics
    myRTVF = RTVF(theConeMosaic, ...
                 theOpticsParams, ...
                 theTargetVisualRFDoGparams, ...
                 'visualizeSpectrallyWeightedPSFs', ~true);
end

function theConeMosaic = generateTestConeMosaic(positionDegs)

    % Set cone aperture modifiers
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
        'smoothLocalVariations', true, ...
        'sigma',  sigmaGaussian, ...
        'shape', 'Gaussian');

    % One L- and one M-cone only
    theConeMosaic = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...
        'eccentricityDegs', positionDegs, ...
        'sizeDegs', [0.2 0.1]*(1+sqrt(sum(positionDegs.^2,2))), ...
        'whichEye', 'right eye', ...
        'customDegsToMMsConversionFunction', @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x), ...
        'customMMsToDegsConversionFunction', @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x), ...
        'overlappingConeFractionForElimination', 0.5, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', coneApertureModifiers);

    

end

