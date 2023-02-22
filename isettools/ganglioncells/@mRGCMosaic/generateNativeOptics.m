function generateNativeOptics(obj, opticsParams)

    % Generate the native optics.
    % These optics form the basis on which surround cone weights are optimized
    % for wiring so as to generate RF with the target visual properties
       
    if (isempty(opticsParams.positionDegs))
        ZernikeCoefficientsEccDegs = obj.eccentricityDegs;
    else
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


    oiEnsemble = obj.inputConeMosaic.oiEnsembleGenerate(...
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

    % Save the native optics params and the native optics
    opticsParams.ZernikeCoefficientsEccDegs = ZernikeCoefficientsEccDegs;
    opticsParams.subtractCentralRefraction = subtractCentralRefraction;
    opticsParams.testSubjectID = testSubjectID;

    obj.theNativeOpticsParams = opticsParams;
    obj.theNativeOptics = oiEnsemble{1};
end


    