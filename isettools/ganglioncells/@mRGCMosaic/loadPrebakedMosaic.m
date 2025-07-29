function [theMRGCMosaic, theOI] = loadPrebakedMosaic(mosaicParams, opticsParams)

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
    
    mRGCMosaicFilename = sprintf('MRGCMosaic_RE_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_Phi_%1.2f_%s_srndModel_%s.mat', ...
        mosaicParams.eccDegs(1), mosaicParams.eccDegs(2), ...
        mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2), ...
        mosaicParams.spatialCompactnessSpectralPurityTradeoff, ...
        opticsSubString, ...
        mosaicParams.surroundOptimizationSubString);

    theFileName = fullfile(prebakedMRGCMosaicDir,mRGCMosaicFilename);
    assert(isfile(theFileName), 'Could not locate the mosaic. File %s not found.\n', theFileName);

    % Load the mosaic 
    load(fullfile(prebakedMRGCMosaicDir,mRGCMosaicFilename), 'theMRGCMosaic')


    if (isfield(mosaicParams, 'cropParams'))&&(~isempty(mosaicParams.cropParams))
        % Crop the mosaic to requested size
        theMRGCMosaic.cropToSizeAtEccentricity(mosaicParams.cropParams.sizeDegs, mosaicParams.cropParams.eccentricityDegs);
    end
    
    % Generate the optics for the mosaic
    [theOI, thePSF] = RGCMosaicAnalyzer.compute.opticsForResponses(...
        theMRGCMosaic, opticsParams.type, ...
        opticsParams.refractiveErrorDiopters, ...
        opticsParams.visualizePSFonTopOfConeMosaic);
 end

