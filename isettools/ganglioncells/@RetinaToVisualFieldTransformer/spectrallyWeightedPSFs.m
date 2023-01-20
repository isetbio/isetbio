function spectrallyWeightedPSFs(obj)

    [obj.theSpectrallyWeightedPSFData, ...
     obj.testSubjectID, ...
     obj.subtractCentralRefraction, ...
     obj.opticsParams] = computeWeightedPSFs(...
            obj.opticsParams, ...
            obj.theConeMosaic, ...
            obj.psfWavelengthSupport);
end


function [thePSFData, testSubjectID, subtractCentralRefraction, opticsParams] = ...
    computeWeightedPSFs(opticsParams, theConeMosaic, psfWavelengthSupport)

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

    % Compute L/M cone weights for weighting the PSF/OTF
    [theLconeSpectralWeights, theMconeSpectralWeights] = LMconeSpectralWeightings(theConeMosaic, opticsParams.positionDegs);

    % L-cone fundamental weighted PSF
    thePSFData.LconeWeighted = spectrallyWeightedPSF(...
        theLconeSpectralWeights, thePSFData.data, ...
        theConeMosaic.wave, psfWavelengthSupport);

    % M-cone fundamental weighted PSF
    thePSFData.MconeWeighted = spectrallyWeightedPSF(...
        theMconeSpectralWeights, thePSFData.data, ...
        theConeMosaic.wave, psfWavelengthSupport);

    % L+M weighted PSF with weights proportional to the number of L and
    % M-cones in the cone mosaic at hand
    lConesNum = numel(find(theConeMosaic.coneTypes == cMosaic.LCONE_ID));
    mConesNum = numel(find(theConeMosaic.coneTypes == cMosaic.MCONE_ID));
    LconePSFweight = lConesNum / (lConesNum  + mConesNum);
    MconePSFweight = mConesNum / (lConesNum  + mConesNum);

    thePSFData.LMconeWeighted = ...
        LconePSFweight * thePSFData.LconeWeighted + ...
        MconePSFweight * thePSFData.MconeWeighted;
    thePSFData.LMconeWeighted = thePSFData.LMconeWeighted / sum(thePSFData.LMconeWeighted(:));
    
    % Specify support in degs instead of the default arc min
    thePSFData.psfSupportXdegs = thePSFData.supportX/60;
    thePSFData.psfSupportYdegs = thePSFData.supportY/60;

    hFig = figure(444); clf;
    set(hFig, 'Position', [10 10 1300 500], 'Color', [1 1 1]);
    subplot(1,3,1);
    imagesc(thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.LconeWeighted);
    axis 'image';
    set(gca, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'FontSize', 16);
    title('L-cone center PSF')

    subplot(1,3,2);
    imagesc(thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.MconeWeighted);
    axis 'image';
    set(gca, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'FontSize', 16);
    title('M-cone center PSF');

    subplot(1,3,3);
    imagesc(thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.LMconeWeighted);
    axis 'image';
    set(gca, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'FontSize', 16);
    title(sprintf('surround PSF (Lw:%2.2f, Mw%2.2f)', LconePSFweight, MconePSFweight));
    colormap(gray(1024))
    drawnow;

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'data');
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
    thePSFData = rmfield(thePSFData, 'supportX');
    thePSFData = rmfield(thePSFData, 'supportY');
end


function [L,M] = LMconeSpectralWeightings(theConeMosaic, theTargetEccDegs)

    if (theConeMosaic.eccVaryingMacularPigmentDensity)
        % Each cone has a different boost factor, based on its eccentricity
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theTargetEccDegs);
    else
        % Single factor based on the eccentricity of the mosaic
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theConeMosaic.eccentricityDegs);
    end
    boostVector = macularPigmentBoostFactors';
    

    % Quantal efficiencies with the foveal MP density and no lens
    coneQuantalEfficienciesFovealMP = theConeMosaic.qe;

    % Quantal efficiencies with the MP density at the target ecc and no lens
    coneQuantalEfficienciesEccBasedMP = diag(boostVector) * coneQuantalEfficienciesFovealMP;

    % Quantal efficiencies with the MP density at the target ecc and the lens
    theLens = Lens('wave',theConeMosaic.wave);
    lensTransmittance = theLens.transmittance;
    coneQuantalEfficienciesFinal = diag(lensTransmittance)*coneQuantalEfficienciesEccBasedMP;

    % Return the L- and M-cone quantal efficiencies
    L = coneQuantalEfficienciesFinal(:,1);
    M = coneQuantalEfficienciesFinal(:,2);
end

function weightedPSF = spectrallyWeightedPSF(spectralWeights, theFullPSF, ...
    theConeMosaicWavelengthSupport, thePSFWavelengthSupport)

    % Compute spectrally-weighted PSF
    weightedPSF = zeros(size(theFullPSF,1), size(theFullPSF,2));
    for iWave = 1:size(theFullPSF,3)
        if (ismember(theConeMosaicWavelengthSupport(iWave), thePSFWavelengthSupport)) || ...
           (isempty(thePSFWavelengthSupport))
            weightedPSF = weightedPSF + theFullPSF(:,:,iWave) * spectralWeights(iWave);
        else
            error('How can this be?');
        end
    end

    % The weighted PSF with unit-volume
    weightedPSF = weightedPSF/sum(weightedPSF(:));
end

