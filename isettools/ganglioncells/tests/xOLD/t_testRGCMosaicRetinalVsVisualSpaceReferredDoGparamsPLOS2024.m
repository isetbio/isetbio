% Script to contrast DoG model params extracted from retinal-space vs visual-space referred STFs
% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicRetinalVsVisualSpaceReferredDoGparamsPLOS2024
%}


    % Optics for STF responses
    if (retinalSpaceReferredSTFs)
        opticsForSTFresponses = 'adaptiveOptics6MM';
    else
        opticsForSTFresponses = 'nativeOptics';
    end

    intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();
    intermediateDataDir = fullfile(intermediateDataDir, 'SLIM')

    % Generate the visual space-referred summary data file
    subjectID = 2;
    subjectID = 3;
    opticsForSTFresponses = 'nativeOptics';
    summaryDataFileName = sprintf('MRGCMosaic_RE_Phi_1.00_Optics_Polans2015-%d_maxStrehlRatio_srndModel_PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0_mRGCMosaic_%s_Achromatic_AllEccentricities.mat', subjectID, opticsForSTFresponses);
    theVisualSpaceReferredSummaryDataFileName = fullfile(intermediateDataDir, 'demos', 'CronerKaplanAnalyses', summaryDataFileName);

    % Generate the retinal space-referred summary data file
    opticsForSTFresponses = 'adaptiveOptics6MM';
    summaryDataFileName = sprintf('MRGCMosaic_RE_Phi_1.00_Optics_Polans2015-%d_maxStrehlRatio_srndModel_PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0_mRGCMosaic_%s_Achromatic_AllEccentricities.mat', subjectID, opticsForSTFresponses);
    theRetinalSpaceReferredSummaryDataFileName = fullfile(intermediateDataDir, 'demos', 'CronerKaplanAnalyses', summaryDataFileName);

    theVisualSpaceReferredDataOut = load(theVisualSpaceReferredSummaryDataFileName)
    theRetinalSpaceReferredDataOut = load(theRetinalSpaceReferredSummaryDataFileName)

    figNo = 1;
    pdfFileName = 'RsRcRatios.pdf';
    RGCMosaicAnalyzer.visualize.correspondenceScatterPlot(figNo, ...
            theVisualSpaceReferredDataOut.RsToRcMosaic, ...
            theRetinalSpaceReferredDataOut.RsToRcMosaic, ...
            2, 'o', [0.5 0.5 0.5], 0.5, 0, [0 20], 0:2:20, 'Rs/Rc', pdfFileName);


    figNo = 2;
    pdfFileName = 'intSCRatios.pdf';
    RGCMosaicAnalyzer.visualize.correspondenceScatterPlot(figNo, ...
            theVisualSpaceReferredDataOut.intStoCsensMosaic, ...
            theRetinalSpaceReferredDataOut.intStoCsensMosaic, ...
            2, 'o', [0.5 0.5 0.5], 0.5, 0, [0 2], 0:0.2:2, 'ISs/ISc', pdfFileName);

