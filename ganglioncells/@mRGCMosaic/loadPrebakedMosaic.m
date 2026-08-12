function [theMRGCMosaic, theOI, thePSF, prebakedMRGCMosaicDir, mRGCMosaicFilename, opticsParams] = loadPrebakedMosaic(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, varargin)
% Load a pre-baked ON-mRGC circuit and, optionally, its optics
%
% Syntax
%   [theMRGCMosaic, theOI, thePSF, prebakedMRGCMosaicDir, ...
%    mRGCMosaicFilename, opticsParams] = mRGCMosaic.loadPrebakedMosaic(...
%        coneMosaicSpecies, opticsSubjectName, rgcMosaicName, ...
%        targetVisualSTFdescriptor, ...)
%
% Returns
%   theMRGCMosaic         - The loaded circuit, empty when none matches
%   theOI, thePSF         - Optics for the circuit, empty when not computed
%   prebakedMRGCMosaicDir - Directory actually holding the loaded file. This
%        is the local gallery for a locally synthesized circuit, and the SDR
%        cache directory for a downloaded one. It does not necessarily join
%        with mRGCMosaicFilename below; use sdrPrebakedCircuitFile to get a
%        full path.
%   mRGCMosaicFilename    - Historical ISETBio filename identifying the
%        circuit. The deposited asset uses a different, normalized name.
%   opticsParams          - Optics parameters for the circuit
%
% Description
%   The circuit is resolved by sdrPrebakedCircuitFile, which prefers a local
%   copy and otherwise downloads the published asset from the isetbio-mosaics
%   SDR deposit and verifies it against the manifest checksum. Circuits are
%   large, so a first load can take some time.
%
% See also
%   sdrPrebakedCircuitFile, mRGCMosaic.listPrebakedMosaics,
%   sdrLegacyFilenameMap

    % Parse input
    p = inputParser;
    p.addParameter('computeTheMosaicOptics', true, @islogical);
    p.addParameter('opticsToEmploy', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('wavefrontSpatialSamples', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('cropParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('onlyReturnMosaicFilename', false, @islogical);
    p.parse(varargin{:});
    
    cropParams = p.Results.cropParams;
    opticsToEmploy = p.Results.opticsToEmploy;
    computeTheMosaicOptics = p.Results.computeTheMosaicOptics;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    onlyReturnMosaicFilename = p.Results.onlyReturnMosaicFilename;
    
    % Generate pStruct with synthesized mosaic params
    pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
        coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor);

    % Extract mosaic and optics params for the original prebaked mosaic
    [mosaicParams, opticsParams] = RGCMosaicConstructor.helper.utils.extractSynthesizedMosaicAndOpticsParams(...
        pStruct);

    % Update mosaicParams with optional cropParams
    if (~isempty(cropParams))
        mosaicParams.cropParams = cropParams;
    end

    % Update opticsParams with optional opticsToEmploy
    if (~isempty(opticsToEmploy))
            opticsParams.type = opticsToEmploy.type;
            if (isfield(opticsToEmploy, 'refractiveErrorDiopters'))
                opticsParams.refractiveErrorDiopters = opticsToEmploy.refractiveErrorDiopters;
            end
    end

    
    switch (opticsParams.ZernikeDataBase)
        case 'Polans2015'
            subjectIndex = find(PolansOptics.constants.subjectRanking == opticsParams.subjectRankOrder);
            opticsSubjectSubstring = sprintf('%s-%d', opticsParams.ZernikeDataBase, subjectIndex);
        case 'Artal2012'
            subjectIndex = find(ArtalOptics.constants.subjectRanking == opticsParams.subjectRankOrder);
            opticsSubjectSubstring = sprintf('%s-%d', opticsParams.ZernikeDataBase, subjectIndex);
        case 'Thibos'
            subjectIndex = find(ThibosOptics.constants.subjectRanking == opticsParams.subjectRankOrder);
            opticsSubjectSubstring = sprintf('%s-%d', opticsParams.ZernikeDataBase, subjectIndex);
        otherwise
            error('Unknown optics database: ''%s''.\n', opticsParams.ZernikeDataBase);
    end % switch

    [~, prebakedMRGCMosaicDir] = mRGCMosaic.listPrebakedMosaics();
    opticsSubString = sprintf('Optics_%s_maxStrehlRatio', opticsSubjectSubstring);
    
    if (isempty(pStruct.customLMSconeDensities))
        mRGCMosaicFilename = sprintf('MRGCMosaic_RE_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_Phi_%1.2f_%s_srndModel_%s.mat', ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2), ...
            mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2), ...
            mosaicParams.spatialCompactnessSpectralPurityTradeoff, ...
            opticsSubString, ...
            mosaicParams.surroundOptimizationSubString);
    else
        mRGCMosaicFilename = sprintf('MRGCMosaic_%1.0f_%1.0f_%1.0f_RE_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_Phi_%1.2f_%s_srndModel_%s.mat', ...
            pStruct.customLMSconeDensities(1)*100, pStruct.customLMSconeDensities(2)*100, pStruct.customLMSconeDensities(3)*100, ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2), ...
            mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2), ...
            mosaicParams.spatialCompactnessSpectralPurityTradeoff, ...
            opticsSubString, ...
            mosaicParams.surroundOptimizationSubString);
    end


    if (onlyReturnMosaicFilename)
        theMRGCMosaic = [];
        theOI = [];
        thePSF = [];
        return;
    end

    % Resolve the circuit locally, or download and verify it from the
    % published isetbio-mosaics SDR deposit.
    theFileName = sdrPrebakedCircuitFile(mRGCMosaicFilename);
    if (isempty(theFileName))
        fprintf(2, 'Could not locate a prebaked mosaic with the provided specifiers.\n');
        fprintf(2, 'Wanted: %s\n', mRGCMosaicFilename);
        fprintf(2, 'Use mRGCMosaic.listPrebakedMosaics to see what is available.\n');
        fprintf(2, 'See ''t_visualizePrebakedMosaicAndOptics'' for a different way to load prebaked mosaics\n');
        theMRGCMosaic = [];
        theOI = [];
        thePSF = [];
        return;
    end
    prebakedMRGCMosaicDir = fileparts(theFileName);

    fprintf('Loading prebaked mRGCmosaic from:\n\t%s\n', mRGCMosaicFilename);

    % Load the mosaic
    load(theFileName, 'theMRGCMosaic')


    fprintf('\t---> The full prebaked mRGC mosaic contains %d mRGCs\n', theMRGCMosaic.rgcsNum);
    if (isfield(mosaicParams, 'cropParams'))&&(~isempty(mosaicParams.cropParams))
   
        if (isfield(mosaicParams.cropParams, 'visualizeSpatialRelationshipToSourceMosaic'))
            visualizeSpatialRelationshipToSourceMosaic = mosaicParams.cropParams.visualizeSpatialRelationshipToSourceMosaic;
        else
            visualizeSpatialRelationshipToSourceMosaic = false;
        end

        % Crop the mosaic to requested size
        theMRGCMosaic.cropToSizeAtEccentricity(...
                mosaicParams.cropParams.sizeDegs, ...
                mosaicParams.cropParams.eccentricityDegs, ...
                'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);
       
    end
    
    fprintf('\t---> The cropped mRGC mosaic contains %d mRGCs\n', theMRGCMosaic.rgcsNum);

    % Generate the optics for the mosaic
    if (computeTheMosaicOptics)

        [theOI, thePSF] = RGCMosaicAnalyzer.compute.opticsForResponses(...
            theMRGCMosaic, ...
            opticsParams.type, ...
            opticsParams.refractiveErrorDiopters, ...
            opticsParams.visualizePSFonTopOfConeMosaic, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples);
    else
        theOI = [];
        thePSF = [];
    end

 end

