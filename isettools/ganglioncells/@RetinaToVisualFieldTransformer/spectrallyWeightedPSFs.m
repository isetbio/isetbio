function spectrallyWeightedPSFs(obj)

    [obj.theSpectrallyWeightedPSFData, ...
     obj.testSubjectID, ...
     obj.subtractCentralRefraction, ...
     obj.opticsParams] = computeWeightedPSFs(...
            obj.opticsParams, ...
            obj.theConeMosaic, ...
            obj.psfWavelengthSupport);
end


function [thePSFData, testSubjectID, subtractCentralRefraction, ...
          opticsParams] = computeWeightedPSFs(...
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
    theVlambda2DegWeights = vLambda2DegsWeights(theConeMosaic.wave);

    % Compute L/M cone weights for weighting the PSF/OTF
    [theLconeWeights, theMconeWeights] = LMconeFundamentals(theConeMosaic.wave);

    % V-lambda weighted PSF
    thePSFData.vLambda2DegWeightedData = spectrallyWeightedPSF(...
        theVlambda2DegWeights, thePSFData.data, ...
        theConeMosaic.wave, psfWavelengthSupport);

    % L-cone fundamental weighted PSF
    thePSFData.LconeWeightedData = spectrallyWeightedPSF(...
        theLconeWeights, thePSFData.data, ...
        theConeMosaic.wave, psfWavelengthSupport);

    % M-cone fundamental weighted PSF
    thePSFData.MconeWeightedData = spectrallyWeightedPSF(...
        theMconeWeights, thePSFData.data, ...
        theConeMosaic.wave, psfWavelengthSupport);

    % Specify support in degs instead of the default arc min
    thePSFData.psfSupportXdegs = thePSFData.supportX/60;
    thePSFData.psfSupportYdegs = thePSFData.supportY/60;

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'data');
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
    thePSFData = rmfield(thePSFData, 'supportX');
    thePSFData = rmfield(thePSFData, 'supportY');
end

function vLambda = vLambdaCIE31Weights(wavelengthSupport)
    load T_xyz1931;
    S = WlsToS(wavelengthSupport(:));
    vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
end

function vLambda = vLambda2DegsWeights(wavelengthSupport)
    load T_CIE_Y2;
    S = WlsToS(wavelengthSupport(:));
    vLambda = SplineCmf(S_CIE_Y2,T_CIE_Y2,S);
end

function [L,M] = LMconeFundamentals(wavelengthSupport)
    load T_cones_ss2
    S = WlsToS(wavelengthSupport(:));
    coneFundamentals = SplineCmf(S_cones_ss2,T_cones_ss2 ,S);
    L = coneFundamentals(1,:);
    M = coneFundamentals(2,:);
end

function weightedPSF = spectrallyWeightedPSF(weights, theFullPSF, ...
    theConeMosaicWavelengthSupport, thePSFWavelengthSupport)

    % Compute spectrally-weighted PSF
    weightedPSF = zeros(size(theFullPSF,1), size(theFullPSF,2));
    for iWave = 1:size(theFullPSF,3)
        if (ismember(theConeMosaicWavelengthSupport(iWave), thePSFWavelengthSupport)) || ...
           (isempty(thePSFWavelengthSupport))
            weightedPSF = weightedPSF + theFullPSF(:,:,iWave) * weights(iWave);
        else
            error('How can this be?');
        end
    end

    % The weighted PSF with unit-volume
    weightedPSF = weightedPSF/sum(weightedPSF(:));
end

