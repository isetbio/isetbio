function performComputeVisuallyProjectedAnatomicalConeRcs(mosaicParams, varargin)


    % Ask the user about what optics to use for computing the input cone
    % mosaic STF responses
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


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

    coneRcDegs = mean(anatomicalConeRcDegs);
    anatomicalConeApertureMap = exp(-(Rdegs/coneRcDegs).^2);
    visuallyProjectedConeApertureMap = conv2(thePSFData.data, anatomicalConeApertureMap, 'same');

    theFittedGaussianEllipsoid = MosaicPoolingOptimizer.fitGaussianEllipsoid(...
        thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, visuallyProjectedConeApertureMap, ...
        varargin);


    size(anatomicalConeApertureMap)
    size(thePSFData.data)

    apertureDataStruct = struct(...
       'positionDegs', [0 0], ...
       'RcDegs', coneRcDegs);

    ff = MSreadyPlot.figureFormat('1x2 large');
    hFig = figure(1); clf;
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    psfRangeDegs = 5/60;
    MSreadyPlot.render2DPSF(theAxes{1,1}, ...
            thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
            thePSFData.data, psfRangeDegs, 'PSF and anatomical cone aperture', ff, ...
            'psfAlpha', 0.1, ...
            'withConeApertureData', apertureDataStruct, ...
            'noXLabel', false, ...
            'noYLabel', false);
    
    imagesc(theAxes{1,2}, thePSFData.psfSupportXdegs*60, thePSFData.psfSupportYdegs*60, visuallyProjectedConeApertureMap);
    axis(theAxes{1,2}, 'image');
    axis(theAxes{1,2}, 'xy');
    set(theAxes{1,2}, 'XLim', psfRangeDegs*60*1.05*[-1 1], 'YLim', psfRangeDegs*60*1.05*[-1 1], 'XTick', -10:1:10, 'YTick', -10:1:10);
    colormap(theAxes{1,2},brewermap(1024, 'greys'));
    title(theAxes{1,2},'visually projected cone aperture');

   pause

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

