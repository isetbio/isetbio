% Phase 10:Covisualize RF center with dendritic RF size
function runPhase10(runParams)

    % Load data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios', 'midgetRGCconnectionMatrix');
     
 
    subregion.center = round(1000*WatsonRGCModel.rhoDegsToMMs(runParams.patchEccDegs));
    subregion.size = round(1000*WatsonRGCModel.rhoDegsToMMs(runParams.patchSizeDegs));     

    [semiAxes,rfCenters] = computeRFsizes(runParams.zLevels, runParams.whichLevelsToContour, ...
            midgetRGCconnectionMatrix, conePositionsMicrons, ...
            RGCRFPositionsMicrons, coneSpacingsMicrons, coneTypes, subregion);
  
    plotlabOBJ = plotlabco();
    plotlabOBJ.applyRecipe(...
        'renderer', 'painters', ...
        'axesBox', 'off', ...
        'colorOrder', [0 0 1; 0.2 0.6 0.5], ...
        'axesTickLength', [0.015 0.01],...
        'axesFontSize', 16, ...
        'figureWidthInches', 12/2, ...
        'figureHeightInches', 9/2);
    
    hFig = figure(222); clf;
    midgetDenditicTreeDiameterArcMinDacey = DaceyData(); hold on;
    scatter(midgetDenditicTreeDiameterArcMinDacey(:,1), midgetDenditicTreeDiameterArcMinDacey(:,2), 100);
    
    set(gca, 'XLim', [0.1 40], 'YLim', [0.5 30], 'YScale', 'log', 'XScale', 'log');
    set(gca, 'XTick', [0.1 0.3 1 3 10 30], 'YTick', [1 3 10 30 100]);
    legend({ ...
            sprintf('Dacey & Petersen (1992), N = %d', size(midgetDenditicTreeDiameterArcMinDacey,1)), ...
            });
    ylabel('\it mRGC dendrite diam. (arc min)');
    xlabel('\it eccentricity (degs)');
    drawnow;
    fName = sprintf('SemiAxes1');
    plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(rootDir, 'exports'));
    
    
    plotlabOBJ.applyRecipe(...
        'renderer', 'painters', ...
        'axesBox', 'off', ...
        'colorOrder', [.9 0.1 0.1; 0 0 1; 0.2 0.6 0.5], ...
        'axesTickLength', [0.015 0.01],...
        'axesFontSize', 22, ...
        'figureWidthInches', 12, ...
        'figureHeightInches', 9);
    
    hFig = figure(222); clf;
    eccNeuronsMicrons = sqrt(sum(rfCenters.^2,2));
    eccNeuronsDegs = WatsonRGCModel.rhoMMsToDegs(eccNeuronsMicrons/1000);
    rfCenterRadiusDegs = WatsonRGCModel.rhoMMsToDegs(mean(semiAxes,2)/1000);
    
    scatter(eccNeuronsDegs, 2*rfCenterRadiusDegs*60, 16); hold on;
    
    midgetDenditicTreeDiameterArcMinDacey = DaceyData();
    scatter(midgetDenditicTreeDiameterArcMinDacey(:,1), midgetDenditicTreeDiameterArcMinDacey(:,2), 100);
    
    midgetRGCDendriticFieldRadiusCowey = CoweyPerryData();
    scatter(midgetRGCDendriticFieldRadiusCowey(:,1), 2*midgetRGCDendriticFieldRadiusCowey(:,2)*60, 144);
    
    set(gca, 'XLim', [0.1 40], 'YLim', [0.5 30], 'YScale', 'log', 'XScale', 'log');
    set(gca, 'XTick', [0.1 0.3 1 3 10 30], 'YTick', [1 3 10 30 100]);
    legend({sprintf('ISETBio (2020), N = %d',numel(eccNeuronsDegs)), ...
            sprintf('Dacey & Petersen (1992), N = %d', size(midgetDenditicTreeDiameterArcMinDacey,1)), ...
            sprintf('Perry & Cowey (1984), N = %d', size(midgetRGCDendriticFieldRadiusCowey,1))});
    ylabel('mRGC dendritic field diameter (arc min)');
    xlabel('eccentricity (degs)');
    drawnow;

    fName = sprintf('SemiAxes2');
    plotlabOBJ.exportFig(hFig, 'png', fName, fullfile(rootDir, 'exports'));
end