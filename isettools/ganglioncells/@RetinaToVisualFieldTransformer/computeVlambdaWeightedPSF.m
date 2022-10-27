function [thePSFData, testSubjectID, subtractCentralRefraction, opticsParams] = computeVlambdaWeightedPSF(...
    opticsParams, theConeMosaic, psfWavelengthSupport)


    % Ensure we have a valid eye specification
    assert(ismember(opticsParams.analyzedEye, {'left eye','right eye'}), ...
        'Invalid analyzed eye specification: ''%s''.', opticsParams.analyzedEye);

    assert(ismember(opticsParams.subjectRankingEye, {'left eye','right eye'}), ...
        'Invalid subject rank eye specification: ''%s''.', opticsParams.subjectRankingEye);


    switch (opticsParams.ZernikeDataBase)
        case RetinaToVisualFieldTransformer.Artal
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(opticsParams.subjectRankingEye);
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                opticsParams.analyzedEye, testSubjectID);

        case RetinaToVisualFieldTransformer.Polans
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

    opticsParams.testSubjectID = testSubjectID;
    opticsParams.zeroCenterPSF = true;
    

    if (opticsParams.refractiveErrorDiopters == -999)
        opticsParams.subtractCentralRefraction = false;
        opticsParams.refractiveErrorDiopters = 0;
    else
        opticsParams.subtractCentralRefraction = subtractCentralRefraction;
    end

    [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);

    if (isempty(oiEnsemble))
        fprintf(2,'Could not generate optics at this eccentricity');
    end

    % Extract the OTF & the PSF
    thePSFData = psfEnsemble{1};

    % Compute v_lambda weights for weigthing the PSF/OTF
    weights = vLambdaWeights(theConeMosaic.wave);

    % Compute vLambda weighted PSF
    vLambdaWeightedPSF = zeros(size(thePSFData.data,1), size(thePSFData.data,2));
    for iWave = 1:size(thePSFData.data,3)
        if (ismember(theConeMosaic.wave(iWave), psfWavelengthSupport)) || ...
           (isempty(psfWavelengthSupport))
            vLambdaWeightedPSF = vLambdaWeightedPSF + thePSFData.data(:,:,iWave) * weights(iWave);
        end
    end
    thePSFData.data = vLambdaWeightedPSF;

    % Center and rotate the PSF
    thePSFData.data = RetinaToVisualFieldTransformer.centerAndRotatePSF(thePSFData.data);

    % Ensure we have a unit volume
    thePSFData.data = thePSFData.data / sum(thePSFData.data(:));

    % Specify support in degs instead of the default arc min
    thePSFData.psfSupportXdegs = thePSFData.supportX/60;
    thePSFData.psfSupportYdegs = thePSFData.supportY/60;

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
    thePSFData = rmfield(thePSFData, 'supportX');
    thePSFData = rmfield(thePSFData, 'supportY');

    % Now generate the circular PSF
%     theCircularPSFData = thePSFData;
%     theCircularPSFData.data = RetinaToVisualFieldTransformer.circularlySymmetricPSF(...
%         thePSFData.data, obj.psfCircularSymmetryMode);
%     obj.theCircularPSFData = theCircularPSFData;

end

function w = vLambdaWeights(wavelengthSupport)
    load T_xyz1931;
    S = WlsToS(wavelengthSupport(:));
    T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
    w = T_vLambda/max(T_vLambda(:));
    w = w / sum(w(:));
end