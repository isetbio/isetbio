function dataOut = generateOptics(obj, opticsParams)

    % Generate optics from the passed params.
    % These optics form the basis on which surround cone weights are optimized
    % for wiring so as to generate RF with the target visual properties
       
    if (isempty(opticsParams.positionDegs))
        % 'Native optics', at mosaic's eccentricity
        ZernikeCoefficientsEccDegs = obj.eccentricityDegs;
    else
        % 'Custom optics', at user-supplied position
        ZernikeCoefficientsEccDegs = opticsParams.positionDegs;
    end

    switch (opticsParams.ZernikeDataBase)
        case ArtalOptics.constants.ZernikeDataBaseName
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(opticsParams.subjectRankingEye);

            if (opticsParams.examinedSubjectRankOrder == 0)
                % Subject 0 results in all-zero Zcoeffs, i.e., AO optics
                testSubjectID = 0;
                subtractCentralRefraction = false;
            else
                testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
                subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    opticsParams.analyzedEye, testSubjectID);
            end


        case PolansOptics.constants.ZernikeDataBaseName
            if (~strcmp(opticsParams.subjectRankingEye, 'right eye'))
                error('Polans measurements exist only for the right eye.');
            end
            rankedSujectIDs = PolansOptics.constants.subjectRanking();
            if (opticsParams.examinedSubjectRankOrder == 0)
                % Subject 0 results in all-zero Zcoeffs, i.e., AO optics
                testSubjectID = 0;
                subtractCentralRefraction = false;
            else
                testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
                subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    testSubjectID);
            end

        otherwise
            error('Unknown zernike database: ''%ss'.', opticsParams.ZernikeDataBase);
    end


    [oiEnsemble, psfEnsemble] = obj.inputConeMosaic.oiEnsembleGenerate(...
        ZernikeCoefficientsEccDegs, ...
       'zernikeDataBase', opticsParams.ZernikeDataBase, ...
       'subjectID', testSubjectID, ...
       'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
       'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
       'noLCA', opticsParams.noLCA, ...
       'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
       'subtractCentralRefraction', subtractCentralRefraction, ...
       'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
       'upsampleFactor', opticsParams.psfUpsampleFactor, ...
       'warningInsteadOfErrorForBadZernikeCoeffs', true);

    % Save the optics params and the generated optics
    opticsParams.ZernikeCoefficientsEccDegs = ZernikeCoefficientsEccDegs;
    opticsParams.subtractCentralRefraction = subtractCentralRefraction;
    opticsParams.testSubjectID = testSubjectID;

    dataOut = struct();
    dataOut.opticsParams = opticsParams;
    dataOut.theOptics = oiEnsemble{1};
    dataOut.thePSFData = psfEnsemble{1};

end

