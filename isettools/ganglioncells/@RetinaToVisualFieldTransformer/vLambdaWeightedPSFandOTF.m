function vLambdaWeightedPSFandOTF(obj)
    
    % Ensure we have a valid eye specification
    assert(ismember(obj.opticsParams.analyzedEye, {'left eye','right eye'}), ...
        'Invalid analyzed eye specification: ''%s''.', obj.opticsParams.analyzedEye);

    assert(ismember(obj.opticsParams.subjectRankingEye, {'left eye','right eye'}), ...
        'Invalid subject rank eye specification: ''%s''.', obj.opticsParams.subjectRankingEye);


    switch (obj.opticsParams.ZernikeDataBase)
        case RetinaToVisualFieldTransformer.Artal
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(obj.opticsParams.subjectRankingEye);
            obj.testSubjectID = rankedSujectIDs(obj.opticsParams.examinedSubjectRankOrder);
            obj.subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                obj.opticsParams.analyzedEye, obj.testSubjectID);

        case RetinaToVisualFieldTransformer.Polans
            if (~strcmp(obj.opticsParams.subjectRankingEye, 'right eye'))
                error('Polans measurements exist only for the right eye.');
            end
            rankedSujectIDs = PolansOptics.constants.subjectRanking();
            obj.testSubjectID = rankedSujectIDs(obj.opticsParams.examinedSubjectRankOrder);
            obj.subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                obj.opticsParams.analyzedEye, obj.testSubjectID);

        otherwise
            error('Unknown zernike database: ''%ss'.', obj.opticsParams.ZernikeDataBase);
    end

    % Compute optics at the rf position
    [oiEnsemble, psfEnsemble] = obj.theConeMosaic.oiEnsembleGenerate(obj.opticsParams.rfPositionEccDegs, ...
                    'zernikeDataBase', obj.opticsParams.ZernikeDataBase, ...
                    'subjectID', obj.testSubjectID, ...
                    'pupilDiameterMM', obj.opticsParams.pupilDiameterMM, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', obj.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', obj.opticsParams.wavefrontSpatialSamples, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);

    if (isempty(oiEnsemble))
        fprintf(2,'Could not generate optics at this eccentricity');
    end

    % Extract the OTF & the PSF
    thePSFData = psfEnsemble{1};
    theOI = oiEnsemble{1};
    theOTFData.data = theOI.optics.OTF.OTF;
    theOTFData.supportX = theOI.optics.OTF.fx;
    theOTFData.supportY = theOI.optics.OTF.fy;
    theOTFData.supportWavelength = theOI.optics.OTF.wave;

    % Compute v_lambda weights for weigthing the PSF/OTF
    weights = vLambdaWeights(obj.theConeMosaic.wave);

    % Compute vLambda weighted PSF
    vLambdaWeightedPSF = zeros(size(thePSFData.data,1), size(thePSFData.data,2));
    for iWave = 1:size(theOTFData.data,3)
        if (ismember(obj.theConeMosaic.wave(iWave), obj.psfWavelengthSupport)) || ...
           (isempty(obj.psfWavelengthSupport))
            vLambdaWeightedPSF = vLambdaWeightedPSF + thePSFData.data(:,:,iWave) * weights(iWave);
        end
    end
    thePSFData.data = vLambdaWeightedPSF;
    thePSFData.data = thePSFData.data / sum(thePSFData.data(:));

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');

    % Now generate the circular PSF
    theCircularPSFData = thePSFData;
    theCircularPSFData.data = RetinaToVisualFieldTransformer.circularlySymmetricPSF(...
        thePSFData.data, obj.psfCircularSymmetryMode);

    obj.thePSFData = thePSFData;
    obj.theCircularPSFData = theCircularPSFData;
end

function w = vLambdaWeights(wavelengthSupport)
    load T_xyz1931;
    S = WlsToS(wavelengthSupport(:));
    T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
    w = T_vLambda/max(T_vLambda(:));
    w = w / sum(w(:));
end