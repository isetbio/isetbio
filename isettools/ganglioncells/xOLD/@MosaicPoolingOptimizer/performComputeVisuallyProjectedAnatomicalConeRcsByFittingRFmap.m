function performComputeVisuallyProjectedAnatomicalConeRcsByFittingRFmap(mosaicParams, varargin)

    retinalMeridian = 'temporal';
    if (isempty(mosaicParams))
        switch (retinalMeridian)
            case 'nasal'
                % Optic disk 13-17
                xEcc = [0.1  0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.00 2.5 3.00 3.50 4.0 4.5 5.0  6.0  7.0  8.0  9.0 10.0 12 19 20 22 24 26 28 30 ];
                yEcc = xEcc*0;
                eccDegs = [xEcc(:) yEcc(:)];
                sizeDegs = [0.05 0.07 0.08 0.09 0.1 0.12 0.13 0.15 0.17 0.2 0.26 0.29 0.3 0.3 0.35 0.35 0.4 0.45 0.45 0.5 0.5 0.5 0.5 0.5 0.5 0.6 0.65 0.65];
                temporalEquivalentEccDegs = temporalEquivalentEccentricityForNasalEccentricity(eccDegs);
            case 'temporal'
                xEcc = -[0.1  0.25 0.50 0.75 1.0     1.25 1.50 1.75 2.00 2.5 3.00 3.50 4.0 4.5 5.0  6.0  7.0  8.0  9.0 10.0 12 14 16 18 20 22 24 26 28 30 ];
                yEcc = xEcc*0;
                eccDegs = [xEcc(:) yEcc(:)];
                sizeDegs = [0.05 0.07 0.09 0.09 0.11 0.12 0.14 0.15 0.17 0.2 0.26 0.29 0.3 0.33 0.35 0.35 0.4 0.45 0.5 0.5 0.5 0.5 0.5 0.6  0.6 0.6 0.6 0.7 0.7 0.7];

                temporalEquivalentEccDegs = eccDegs;
            otherwise
                error('meridian %s not coded', retinalMeridian);
        end


%         
        for iEcc = 1:size(eccDegs,1)
            mosaicParams = struct(...
                'eccDegs', eccDegs(iEcc,:), ...
                'sizeDegs', sizeDegs(iEcc)*[1 1], ...
                 'rgcType', 'ONcenterMidgetRGC');

            
            opticsParams = mRGCMosaic.defaultOpticsParams;
            if (abs(mosaicParams.eccDegs(1))> 15)
                opticsParams.wavefrontSpatialSamples = 901;
            end
            
       
            [visuallyProjectedConeAperturesRcDegs, coneRFpositionsDegs] = doIt(mosaicParams, opticsParams);

            theVisuallyProjectedConeApertureRcDegs(iEcc) = median(visuallyProjectedConeAperturesRcDegs);

            figure(10); clf;
            plot(abs(temporalEquivalentEccDegs(1:iEcc,1)), theVisuallyProjectedConeApertureRcDegs(1:iEcc), 'rs-');
            set(gca, 'XLim', [0.05 30], 'XScale', 'log', 'YLim', [0 0.1], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], 'YTick', 0:0.01:0.1);
            drawnow;
        end

        save(sprintf('visuallyProjectedConeApertures%sRetina.mat', retinalMeridian), 'eccDegs', 'temporalEquivalentEccDegs', 'theVisuallyProjectedConeApertureRcDegs');

    else
        % Ask the user about what optics to use for computing the input cone
        % mosaic STF responses
        opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


        [theVisuallyProjectedConeApertureRcDegs, coneRFpositionsDegs] = doIt(mosaicParams, opticsParams);
    end


end

function [theVisuallyProjectedConeApertureRcDegs, coneRFpositionsDegs] = doIt(mosaicParams, opticsParams)

    
    % Generate input cone mosaic
    inputConeMosaic = generateInputConeMosaic(mosaicParams, opticsParams);

    % Compute anatomical cone Rc
    anatomicalConeRcDegs = inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                           inputConeMosaic.coneApertureDiametersDegs;

    % Cone positions
    coneRFpositionsDegs = inputConeMosaic.coneRFpositionsDegs;


    % Optics params
    if (isempty(opticsParams.positionDegs))
        % 'Native optics', at mosaic's eccentricity
        ZernikeCoefficientsEccDegs = mosaicParams.eccDegs;
    else
        % 'Custom optics', at user-supplied position
        ZernikeCoefficientsEccDegs = opticsParams.positionDegs;
    end

    switch (opticsParams.ZernikeDataBase)
        case ArtalOptics.constants.ZernikeDataBaseName
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(opticsParams.subjectRankingEye);
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                opticsParams.analyzedEye, testSubjectID);

        case PolansOptics.constants.ZernikeDataBaseName
            if (~strcmp(opticsParams.subjectRankingEye, 'right eye'))
                error('Polans measurements exist only for the right eye.');
            end
            rankedSujectIDs = PolansOptics.constants.subjectRanking();
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                testSubjectID);

        otherwise
            error('Unknown zernike database: ''%ss'.', opticsParams.ZernikeDataBase);
    end

    % Generate optics
    [oiEnsemble, psfEnsemble] = inputConeMosaic.oiEnsembleGenerate(...
        ZernikeCoefficientsEccDegs, ...
       'zernikeDataBase', opticsParams.ZernikeDataBase, ...
       'subjectID', testSubjectID, ...
       'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
       'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
       'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
       'subtractCentralRefraction', subtractCentralRefraction, ...
       'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
       'upsampleFactor', opticsParams.psfUpsampleFactor, ...
       'warningInsteadOfErrorForBadZernikeCoeffs', true);

    % Retrieve the optics and the PSF
    theOptics = oiEnsemble{1};
    thePSFData = psfEnsemble{1};

    % Transform support from arc min to degs
    thePSFData.psfSupportXdegs = thePSFData.supportX/60;
    thePSFData.psfSupportYdegs = thePSFData.supportY/60;

    targetWavelength = 550;
    [~,idx] = min(abs(inputConeMosaic.wave-targetWavelength));
    thePSFData.data = squeeze(thePSFData.data(:,:,idx));

    % Flip left-right
    thePSFData.data = fliplr(thePSFData.data);

    [Xdegs,Ydegs] = meshgrid(thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs);
    Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

    theVisuallyProjectedConeApertureRcDegs = zeros(1, numel(anatomicalConeRcDegs));

    % Fit each cone
    visualizeFits = true;
    fprintf('Ftting apertures of %d cones\n', numel(anatomicalConeRcDegs));
    for iCone = 1:numel(anatomicalConeRcDegs)
        if (iCone > 5)
            visualizeFits = false;
        end

        fprintf('Ftting aperture of cone %d of %d (ecc: %2.2f,%2.2f)\n', iCone, numel(anatomicalConeRcDegs), coneRFpositionsDegs(iCone,1), coneRFpositionsDegs(iCone,2));
        theVisuallyProjectedConeApertureRcDegs(iCone) = ...
            computeVisuallyProjectedConeAperture(iCone,anatomicalConeRcDegs(iCone), Rdegs, thePSFData,visualizeFits);
    end

end

function theVisuallyProjectedConeApertureRcDegs = computeVisuallyProjectedConeAperture(iCone,coneRcDegs, Rdegs, thePSFData, visualizeFits)

    anatomicalConeApertureMap = exp(-(Rdegs/coneRcDegs).^2);
    visuallyProjectedConeApertureMap = conv2(thePSFData.data, anatomicalConeApertureMap, 'same');

    theFittedGaussian = MosaicPoolingOptimizer.fitGaussianEllipsoid(...
        thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, visuallyProjectedConeApertureMap);

    

    if (visualizeFits)
        ff = MSreadyPlot.figureFormat('1x2 large');
        hFig = figure(iCone); clf;
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);
    
        apertureDataStruct = struct(...
            'positionDegs', [0 0], ...
            'RcDegs', coneRcDegs);

        psfRangeDegs = 5/60;
        MSreadyPlot.render2DPSF(theAxes{1,1}, ...
                thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                thePSFData.data, psfRangeDegs, 'PSF and anatomical cone aperture', ff, ...
                'psfAlpha', 0.1, ...
                'withConeApertureData', apertureDataStruct, ...
                'noXLabel', false, ...
                'noYLabel', false);
        
        imagesc(theAxes{1,2}, thePSFData.psfSupportXdegs*60, thePSFData.psfSupportYdegs*60, visuallyProjectedConeApertureMap);
        hold(theAxes{1,2}, 'on');
        zLevels = max(theFittedGaussian.ellipsoidMap(:))*exp(-1.0)*[0.99 1.0];
        contour(thePSFData.psfSupportXdegs*60, thePSFData.psfSupportYdegs*60, theFittedGaussian.ellipsoidMap, zLevels, ...
            'Color', [1 0 0], 'LineWidth' , 1.5);
    
        axis(theAxes{1,2}, 'image');
        axis(theAxes{1,2}, 'xy');
        set(theAxes{1,2}, 'XLim', psfRangeDegs*60*1.05*[-1 1], 'YLim', psfRangeDegs*60*1.05*[-1 1], 'XTick', -10:1:10, 'YTick', -10:1:10);
        colormap(theAxes{1,2},brewermap(1024, 'greys'));
        title(theAxes{1,2}, sprintf('visually projected cone aperture, Rc (arc min)= %2.2f, %2.2f', ...
            60*min(theFittedGaussian.characteristicRadii), 60*max(theFittedGaussian.characteristicRadii)));

        %theFittedGaussian.ellipsoidMap = gaussian2D(fittedParams,xydata);
        %theFittedGaussian.characteristicRadii = [RcX RcY];
        %theFittedGaussian.rotationDegs = fittedParams(6);
        %theFittedGaussian.flatTopExponents = [fittedParams(7) fittedParams(8)];
        %theFittedGaussian.xyCenter = [xo yo];

    end

    

    theVisuallyProjectedConeApertureRcDegs = min(theFittedGaussian.characteristicRadii);
end






function inputConeMosaic = generateInputConeMosaic(mosaicParams, opticsParams)

   % Set cone aperture modifiers
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
        'smoothLocalVariations', true, ...
        'sigma',  sigmaGaussian, ...
        'shape', 'Gaussian');
    
    % Generate the input cone mosaic
    inputConeMosaic = cMosaic(...
        'sourceLatticeSizeDegs', 58, ...
        'eccentricityDegs', mosaicParams.eccDegs, ...
        'sizeDegs', mosaicParams.sizeDegs, ...
        'whichEye', opticsParams.analyzedEye, ...
        'eccVaryingConeBlur', true, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'overlappingConeFractionForElimination', 0.5, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', coneApertureModifiers);

end


function temporalEquivalentEccDegs = temporalEquivalentEccentricityForNasalEccentricity(eccDegs)

    if (size(eccDegs,2) ~= 2)
        error('eccDegs must be an N x 2 matrix of (x,y) eccentricities.');
    end
    dataPoints = size(eccDegs,1);

    temporalEquivalentEccDegs = 0*eccDegs;
    
    for iPoint = 1:dataPoints
        % The vertical ecc
        verticalEcc = eccDegs(iPoint,2);
        horizontalEcc = 0.70 * eccDegs(iPoint,1);

        % The radial temporal equivalent eccentricity
        radialEquivalentEccentricity = sqrt(horizontalEcc^2 + verticalEcc^2);

        % Compute the vector temporal equivalent eccentricity
        theta = angle(eccDegs(iPoint,1) + 1j*eccDegs(iPoint,2));
        temporalEquivalentEccDegs(iPoint,:) = radialEquivalentEccentricity  * [-cos(theta) sin(theta)];
    end % iPoint
end
