function computeSpectrallyWeightedPSFs(obj, visualize, varargin)
% Compute spectrally weighted (L-cone, M-cone and L+M-cone weighted) PSFs
%
% Syntax:
%   spectrallyWeightedPSFs(obj, visualize, varargin)
%
% Description:
%    Computes spectrally weighted (L-cone, M-cone and L+M-cone weighted)
%    PSFs, where the L- and M-w spectral eights are derived from the retinal
%    spectral absorptance of actual L- and M-cones in the input cone
%    mosaic at the desrired eccentricity (opticsParams.positionDegs), 
%    taking into account macular pigment density at that location, as 
%    specified by the input @cMosaic.
%    This function sets the obj.theSpectrallyWeightedPSFData parameter
%    and updates the obj.opticsParams
%
% Inputs:
%    obj                 - An @RTVF object
%
% Outputs:
%    none
%
% Optional key/value pairs:
%    none
%         

    % Parse input
    %p = inputParser;
    %p.parse(varargin{:});

    [obj.spectrallyWeightedPSFData, ...
     obj.opticsParams] = computeSpecrallyWeightedPSFs(...
            obj.opticsParams, ...
            obj.coneMosaic, ...
            visualize);
end


function [thePSFData, opticsParams] = computeSpecrallyWeightedPSFs(opticsParams, theConeMosaic, visualize)
    % Ensure we have a valid eye specification
    assert(ismember(opticsParams.analyzedEye, {'left eye','right eye'}), ...
        'Invalid analyzed eye specification: ''%s''.', opticsParams.analyzedEye);

    assert(ismember(opticsParams.subjectRankingEye, {'left eye','right eye'}), ...
        'Invalid subject rank eye specification: ''%s''.', opticsParams.subjectRankingEye);


    switch (opticsParams.ZernikeDataBase)
        case RTVF.Artal
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(opticsParams.subjectRankingEye);
            testSubjectID = rankedSujectIDs(opticsParams.examinedSubjectRankOrder);
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                opticsParams.analyzedEye, testSubjectID);

        case RTVF.Polans
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

    eccDegs = sqrt(sum(opticsParams.positionDegs.^2,2));
    if (eccDegs < 0.5) && (opticsParams.psfUpsampleFactor < 2)
        opticsParams.psfUpsampleFactor = 2;
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
    [theLconeSpectralWeights, theMconeSpectralWeights] = ...
        RTVF.LMconeSpectralWeightings(theConeMosaic, opticsParams.positionDegs);

    % L-cone fundamental weighted PSF
    thePSFData.LconeWeighted = spectrallyWeightedPSF(...
        theLconeSpectralWeights, thePSFData.data);

    % M-cone fundamental weighted PSF
    thePSFData.MconeWeighted = spectrallyWeightedPSF(...
        theMconeSpectralWeights, thePSFData.data);

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

    % Spatial surrport for computing RF maps. This will be changed when
    % cropping the PSF
    thePSFData.spatialSupportForRFmapXdegs = thePSFData.psfSupportXdegs;
    thePSFData.spatialSupportForRFmapYdegs = thePSFData.psfSupportYdegs;


    if (visualize)
        ff = MSreadyPlot.figureFormat('1x4');
        hFig = figure(4000); clf;
        set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
        
        psfRange = 0.1;

        
        % Render the L-cone weighted PSF
        ax = subplot('Position', ff.subplotPosVectors(1,1).v);
        MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.LconeWeighted, psfRange, 'L-cone PSF (center)', ff, ...
            'noXLabel', false, 'noYLabel', false);
        
        % Render the M-cone weighted PSF
        ax = subplot('Position', ff.subplotPosVectors(1,2).v);
        MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.MconeWeighted, psfRange, 'M-cone PSF (center)', ff, ...
            'noXLabel', false, 'noYLabel', true);
        

        % Render the L+M-cone weighted PSF
        ax = subplot('Position', ff.subplotPosVectors(1,3).v);
        MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.LMconeWeighted, psfRange, 'L+M-cone PSF (surround)', ff, ...
            'noXLabel', false, 'noYLabel', true);
        
        % Render the spectral weights
        ax = subplot('Position', ff.subplotPosVectors(1,4).v);
        MSreadyPlot.renderConeFundamentals(ax, theConeMosaic.wave, theLconeSpectralWeights, theMconeSpectralWeights, [], 'spectral weights', ff, ...
            'noYLabel', true);
        

        drawnow;
    end

    % Remove irrelevant fields
    thePSFData = rmfield(thePSFData, 'data');
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
    thePSFData = rmfield(thePSFData, 'supportX');
    thePSFData = rmfield(thePSFData, 'supportY');
end




function weightedPSF = spectrallyWeightedPSF(spectralWeights, theFullPSF)

    % Compute spectrally-weighted PSF
    weightedPSF = zeros(size(theFullPSF,1), size(theFullPSF,2));
    for iWave = 1:size(theFullPSF,3)
        weightedPSF = weightedPSF + theFullPSF(:,:,iWave) * spectralWeights(iWave);
    end

    % The weighted PSF with unit-volume
    weightedPSF = weightedPSF/sum(weightedPSF(:));
end



