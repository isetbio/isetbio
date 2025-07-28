
%% Introduction to the 2025 midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - how to load of of pre-baked midget RGC mosaics,
%


% History:
%    07/28/25  NPC  Wrote it.

function t_mRGCMosaicBasic
    %% Close all figures
    close all;

    % Load an mRGCmosaic located the far periphery
    theMRGCmosaic = farPeripheryMRGCmosaic();
    theMRGCmosaic.visualize();

end

% Supporting functions
function theMRGCMosaic = farPeripheryMRGCmosaic()

    theOpticsSubject = 'Polans2015-2';
    theMosaicXYeccentricityDegs = [-32.0 0.0];
    theMosaicXYsizeDegs = [9 9];

    prebakedMRGCMosaicDir = 'isettools/ganglioncells/data/prebakedRGCmosaics/ONmRGCmosaics';
    spatialCompactnessSpectralPurityTradeoff = 1;
    opticsSubString = sprintf('Optics_%s_maxStrehlRatio', theOpticsSubject);
    surroundOptimizationSubString = 'PackerDacey2002H1freeUpperH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';

    mRGCMosaicFilename = sprintf('MRGCMosaic_RE_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_Phi_%1.2f_%s_srndModel_%s.mat', ...
        theMosaicXYeccentricityDegs(1), theMosaicXYeccentricityDegs(2), ...
        theMosaicXYsizeDegs(1), theMosaicXYsizeDegs(2), ...
        spatialCompactnessSpectralPurityTradeoff, opticsSubString, surroundOptimizationSubString);

    load(fullfile(isetbioRootPath, prebakedMRGCMosaicDir,mRGCMosaicFilename), 'theMRGCMosaic');

    % Employ the native optics (what was used to optimize the surround)
    opticsForSTFresponses = 'nativeOptics';
    %opticsForSTFresponses = 'adaptiveOptics6MM';
    residualWithRespectToNativeOpticsDefocusDiopters = [];
    visualizePSFonTopOfConeMosaic = true;

    % Generate the optics for the mosaic
    [theOI, thePSF] = RGCMosaicAnalyzer.compute.opticsForResponses(...
        theMRGCMosaic, opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, visualizePSFonTopOfConeMosaic);

    
 end


