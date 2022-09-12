function visualizeResults(obj)
    hFig = visualizeTargetAndFittedRFs(obj);
    pdfFileName = strrep(obj.computedObjDataFileName, '.mat', '.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
end

function hFig = visualizeTargetAndFittedRFs(obj)

   

    figNo = 1;
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [15 15 1350 1000], 'Color', [1 1 1]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 4, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.08, ...
           'topMargin',      0.02);

    obj.rfComputeStruct.theFittedVisualRFMap = obj.rfComputeStruct.theFittedVisualRFMap/max(obj.rfComputeStruct.theFittedVisualRFMap(:));
    obj.rfComputeStruct.targetVisualRFMap = obj.rfComputeStruct.targetVisualRFMap / max(obj.rfComputeStruct.targetVisualRFMap(:));
    
    % Extract needed info
    maxPSF = max([max(obj.thePSFData.data(:)) max(obj.theCircularPSFData.data(:))]);
    maxRetinalRF  = max(obj.rfComputeStruct.theRetinalRFsurroundConeMap(:));
    maxRetinalRFprofile = 3*max(sum(obj.rfComputeStruct.theRetinalRFsurroundConeMap,1));

    maxRF = max(obj.rfComputeStruct.targetVisualRFMap(:));
    maxRFprofile = 0.3*max(sum(obj.rfComputeStruct.targetVisualRFMap,1));
    maxVisualRF  = 3*max(obj.rfComputeStruct.targetVisualRFsurroundMap(:));
    maxVisualRFprofile = 3*max(sum(obj.rfComputeStruct.targetVisualRFsurroundMap,1));
    maxFittedRF  = 3*max(obj.rfComputeStruct.theFittedVisualRFsurroundConeMap(:));
    maxFittedRFprofile = 3*max(sum(obj.rfComputeStruct.theFittedVisualRFsurroundConeMap,1));
    conesNumInRFcenter = numel(obj.rfComputeStruct.modelConstants.indicesOfCenterCones);

    maxSpatialSupportDegs = max(obj.rfComputeStruct.modelConstants.spatialSupportDegs(:));
    spatialSupportXdegs = obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1);
    spatialSupportYdegs = obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,2);

    retinalSCratio = sum(obj.rfComputeStruct.theRetinalRFsurroundConeMap(:))/sum(obj.rfComputeStruct.theRetinalRFcenterConeMap(:));
    targetSCratio = sum(obj.rfComputeStruct.targetVisualRFsurroundMap(:))/sum(obj.rfComputeStruct.targetVisualRFcenterMap(:));
    achievedSCratio = sum(obj.rfComputeStruct.theFittedVisualRFsurroundConeMap(:))/sum(obj.rfComputeStruct.theFittedVisualRFcenterConeMap(:));

    eccDegs = obj.theConeMosaic.eccentricityDegs;
    subjectID = obj.testSubjectID;
   

    % The original PSF
    ax = subplot('Position', subplotPosVectors(1,1).v);
    plotPSF(ax, obj.thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, subjectID, true, false, '');

    ax = subplot('Position', subplotPosVectors(1,2).v);
    plotPSF(ax, obj.theCircularPSFData, maxPSF, maxSpatialSupportDegs, eccDegs, subjectID, true, true, 'circularly-averaged PSF');

    % The retinal RF center cone map
    ax = subplot('Position', subplotPosVectors(1,3).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('retinal RF center (%d cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('retinal RF center (%d cones)', conesNumInRFcenter);
    end
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theRetinalRFcenterConeMap, maxRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, true, true);

    % The retinal RF surround cone map
    ax = subplot('Position', subplotPosVectors(1,4).v);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        -obj.rfComputeStruct.theRetinalRFsurroundConeMap, maxRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, sprintf('retinal RF surround, S/Cratio: %2.2f', retinalSCratio), true, true);


    % The target RF center cone map
    ax = subplot('Position', subplotPosVectors(2,1).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('target visual RF center (%d cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('target visual RF center (%d cones)', conesNumInRFcenter);
    end
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.targetVisualRFcenterMap, maxVisualRF, maxVisualRFprofile, ...
        maxSpatialSupportDegs, plotTitle, true,false);

    % The target RF surround cone map
    ax = subplot('Position', subplotPosVectors(2,2).v);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        -obj.rfComputeStruct.targetVisualRFsurroundMap, maxVisualRF, maxVisualRFprofile, ...
        maxSpatialSupportDegs, 'target visual RF surround', true, true);


    % The achieved RF center cone map
    ax = subplot('Position', subplotPosVectors(2,3).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('achieved visual RF center (%d cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('achieved visual RF center (%d cones)', conesNumInRFcenter);
    end
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theFittedVisualRFcenterConeMap, maxFittedRF, maxFittedRFprofile, ...
        maxSpatialSupportDegs, plotTitle, true, true);

    % The achieved RF surround cone map
    ax = subplot('Position', subplotPosVectors(2,4).v);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        -obj.rfComputeStruct.theFittedVisualRFsurroundConeMap, maxFittedRF, maxFittedRFprofile, ...
        maxSpatialSupportDegs, 'achieved visual RF surround', true, true);


    % The target RF cone map
    ax = subplot('Position', subplotPosVectors(3,1).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('target visual RF (%d cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('target visual RF (%d cones)', conesNumInRFcenter);
    end
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.targetVisualRFMap, maxRF, maxRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, false);

    % The achieved RF cone map
    ax = subplot('Position', subplotPosVectors(3,2).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('achieved visual RF (%d cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('achieved visual RF (%d cones)', conesNumInRFcenter);
    end
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theFittedVisualRFMap, maxRF, maxRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, true);

    ax = subplot('Position', subplotPosVectors(3,3).v);
    RMSE = sqrt(1/numel(obj.rfComputeStruct.targetVisualRFMap)*sum(((obj.rfComputeStruct.targetVisualRFMap(:)-obj.rfComputeStruct.theFittedVisualRFMap(:))).^2));
    plotTitle = sprintf('residual, RMSE: %2.3f 1E-3', 1000*RMSE);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.targetVisualRFMap-obj.rfComputeStruct.theFittedVisualRFMap, ...
        maxRF, maxRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, true);

    ax = subplot('Position', subplotPosVectors(3,4).v);
    plotRFprofiles(ax, spatialSupportXdegs, obj.rfComputeStruct.theFittedVisualRFMap, ...
        obj.rfComputeStruct.targetVisualRFMap, maxSpatialSupportDegs, targetSCratio, achievedSCratio);

    % Finally, plot the surround correction factors
    if (strcmp(obj.targetVisualRFDoGparams.retinalConePoolingModel, ...
              'arbitrary center cone weights, gaussian surround weights with adjustments'))
        correctionFactorsSupportDegs = obj.rfComputeStruct.modelConstants.arbitrarySurroundCorrectionRadialSupportDegs;
        correctionFactors = obj.rfComputeStruct.retinalConePoolingParams(4:end);
        plotSurroundCorrectionFactors(correctionFactorsSupportDegs, correctionFactors, maxSpatialSupportDegs);
     end

end

function plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID, noXTickLabel, noYTickLabel, plotTitle)
    psfZLevels = 0.05:0.1:0.95;
    
    contourf(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data/maxPSF, psfZLevels);
    hold on;
    midRow = (size(thePSFData.data,1)-1)/2+1;
    plot(ax, thePSFData.supportXdegs, -maxSpatialSupportDegs + 1.2*thePSFData.data(midRow,:)/maxPSF*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    axis(ax,'image'); axis 'xy';
    
    ticks = ticksForSpatialSupport(maxSpatialSupportDegs);
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
            'XTick', ticks, 'YTick', ticks, 'CLim', [0 1], 'FontSize', 16);
    set(ax, 'XTickLabel', {});

    if (noXTickLabel)
        set(ax, 'XTickLabel', {});
    else
        xlabel('degrees');
    end

    if (noYTickLabel)
        set(ax, 'YTickLabel', {});
    else
        ylabel('degrees');
    end


    grid(ax, 'on');
    xtickangle(ax, 90);
    if (isempty(plotTitle))
        title(ax, sprintf('PSF (%2.2f, %2.2f degs), subj: %d', eccDegs(1), eccDegs(2), testSubjectID));
    else
        title(ax, plotTitle);
    end
    colormap(ax,brewermap(1024, 'greys'));
end


function plotRF(ax, rfSupportX, rfSupportY, RF, maxRF, maxProfile, maxSpatialSupportDegs, ...
    titleString, noXTickLabels, noYTickLabels)

    imagesc(ax, rfSupportX, rfSupportY, RF/maxRF);
    hold on;
    theProfile = sum(RF,1)/maxProfile;
    
    plot(ax, rfSupportX, -maxSpatialSupportDegs*0.75 + 0.5*theProfile*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    plot(ax, rfSupportX, -maxSpatialSupportDegs*0.75 + rfSupportX*0, 'k-');
    axis(ax,'image'); axis 'xy';
    
    ticks = ticksForSpatialSupport(maxSpatialSupportDegs);
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
            'XTick', ticks, 'YTick', ticks, 'CLim', [-1 1], 'FontSize', 16);
    if (noYTickLabels)
        set(ax, 'YTickLabel', {});
    else
        ylabel(ax, 'degrees');
    end

    if (noXTickLabels)
        set(ax, 'XTickLabel', {});
    else
        xlabel(ax,'degrees');
    end


    grid(ax, 'on');
    xtickangle(ax, 90);
    colormap(ax,brewermap(1024, '*RdBu'));
    title(ax, titleString);
end

function plotRFprofiles(ax, rfSupportX, achievedRF, targetRF, maxSpatialSupportDegs, targetSCratio, achievedSCratio)
    
    theProfileX = sum(achievedRF,1);
    theProfileY = sum(achievedRF,2);
    theTargetProfile = sum(targetRF,1);
    
    maxProfile = max([max(abs(theProfileX)) max(abs(theProfileY)) max(abs(theTargetProfile))]);
    theProfileX = theProfileX / maxProfile;
    theProfileY = theProfileY / maxProfile;
    theTargetProfile = theTargetProfile / maxProfile;
    
    shadedAreaPlot(ax,rfSupportX, theTargetProfile, 0, [0.8 0.8 0.8], [0.4 0.4 0.4], 0.7, 1, '-');
    hold(ax, 'on');
    plot(ax, rfSupportX, theProfileX, 'r-', 'LineWidth', 1.5);
    plot(ax, rfSupportX, theProfileY, 'b-', 'LineWidth', 1.5);
    plot(ax, rfSupportX, theTargetProfile-theProfileX, 'k--', 'LineWidth', 1.5);
    
    ticks = ticksForSpatialSupport(maxSpatialSupportDegs);
    axis(ax, 'square');
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', [-0.3 1.01], ...
            'XTick', ticks, 'YTick', -1:0.2:1, 'FontSize', 16);
    
    grid(ax, 'on');
    xtickangle(ax, 90);
    legend({'target', 'X', 'Y', 'residual'});
    title(sprintf('S/C int. ratio - target: %2.2f\nS/C int. ratio - achieved: %2.2f',targetSCratio, achievedSCratio));
    xlabel(ax,'degrees');
end


function plotSurroundCorrectionFactors(correctionFactorsSupportDegs, correctionFactors, maxSpatialSupportDegs)
    figNo = 2;
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [15 15 1350 1000], 'Color', [1 1 1]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 4, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.04, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.02);

    ax = subplot('Position', subplotPosVectors(1,1).v);
        
    % The surround correction factors
    x = [-fliplr(correctionFactorsSupportDegs)  correctionFactorsSupportDegs(2:end)];
    y = [ fliplr(correctionFactors)            correctionFactors(2:end)];
    stem(ax, x, -y, ...
        'LineStyle','-',...
        'MarkerFaceColor',[1 0.5 0.5],...
        'MarkerEdgeColor',[1 0 0]);

    tickSeparationDegs = 0.05;
    ticks = 0:tickSeparationDegs:5;
    ticks = [-fliplr(ticks) ticks(2:end)];

    axis 'square';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', [-0.3 0.3], 'YTick', -1:0.05:1, ...
            'XTick', ticks,  'FontSize', 16);

    grid(ax, 'on');
    xtickangle(ax, 90);
    xlabel(ax,'degrees');
    ylabel(ax,'correction');
    title(ax, 'surround correction factors')
end

function ticks = ticksForSpatialSupport(maxSpatialSupportDegs)
    if (maxSpatialSupportDegs < 0.2)
        tickSeparationDegs = 0.05;
    elseif (maxSpatialSupportDegs < 0.4)
        tickSeparationDegs = 0.1;
    elseif (maxSpatialSupportDegs < 0.6)
        tickSeparationDegs = 0.15;
    else
        tickSeparationDegs = 0.2;
    end

    ticks = 0:tickSeparationDegs:5;
    ticks = [-fliplr(ticks) ticks(2:end)];
end

function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
