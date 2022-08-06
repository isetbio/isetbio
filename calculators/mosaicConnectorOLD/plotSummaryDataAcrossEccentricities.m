function plotSummaryDataAcrossEccentricities
    rgcMosaicPatchHorizontalEccMicrons = [150 300 600 900 1500 3000 4500]; % 6000];
    rgcMosaicPatchSizeMicrons = [30 50 60 70 100 200 250]; % 300];
  
    rootDir = fileparts(which(mfilename()));
    exportsDir = fullfile(rootDir, 'exports');
    exportsDir = fullfile(exportsDir,'LMSConeMosaic');
    
    eccDegs = [];
    LMbalance = [];
    centerRadii = [];
    surroundRadii = [];
    centerSensitivity = [];
    surroundSensitivity = [];
    
    for patchEccIndex = 1:numel(rgcMosaicPatchHorizontalEccMicrons)

        matFile = fullfile(sprintf('%s_HorizontalEcc_%2.0fmicrons',exportsDir, ...
            rgcMosaicPatchHorizontalEccMicrons(patchEccIndex)), ...
            'ResponseDerivedParams_LMS_0.10_0.10_0.00_PolansSID_10_normalOptics.mat');
        
        fprintf('Loading data for %2.1f degs\n', rgcMosaicPatchHorizontalEccMicrons(patchEccIndex));
        
        
        load(matFile, ...
            'eccRadiusDegs', 'LMconeBalance', ...
            'centerCharacteristicRadii', 'surroundCharacteristicRadii', ...
            'centerPeakSensitivities', 'surroundPeakSensitivities');
        eccDegs = cat(1, eccDegs, eccRadiusDegs(:));
        LMbalance = cat(1, LMbalance, LMconeBalance(:));
        centerRadii = cat(1, centerRadii, centerCharacteristicRadii(:));
        surroundRadii = cat(1, surroundRadii, surroundCharacteristicRadii(:));
        centerSensitivity = cat(1,centerSensitivity, centerPeakSensitivities(:));
        surroundSensitivity = cat(1,surroundSensitivity, surroundPeakSensitivities(:));
    end
    
    plotlabOBJ = setupPlotLab(0, 24,15);
    RGCconeInputInfo = [];
    figExportsDir = pwd();
    
    visualizeRFparamsForConnectedPatch(555, 'ResponseDerivedParamsSummary', ...
        RGCconeInputInfo, ...
        eccDegs, ...
        centerRadii, surroundRadii, ...
        centerSensitivity, surroundSensitivity, ...
        sprintf('ResponseDerivedParamsSummary'), ...
        figExportsDir, plotlabOBJ);

    %setupPlotLab(-1);
end

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 

