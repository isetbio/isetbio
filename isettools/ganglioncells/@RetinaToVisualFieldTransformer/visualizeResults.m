function visualizeResults(obj)
    hFig = visualizeTargetAndFittedRFs(obj);
    pdfFileName = strrep(obj.computedObjDataFileName, '.mat', '.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
end

function hFig = visualizeTargetAndFittedRFs(obj)

    % Extract data
    eccDegs = obj.theConeMosaic.eccentricityDegs;
    subjectID = obj.testSubjectID;
    conesNumInRFcenter = numel(obj.rfComputeStruct.modelConstants.indicesOfCenterCones);

    maxSpatialSupportDegs = max(obj.rfComputeStruct.modelConstants.spatialSupportDegs(:));
    spatialSupportXdegs = obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1);
    spatialSupportYdegs = obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,2);
    
    switch (obj.targetVisualRFDoGparams.visualRFmodel)
        case 'gaussian center, gaussian surround'
            visualRFcenterCharacteristicRadiusDegs = obj.rfComputeStruct.targetVisualRFparamsVector(2);

        case 'ellipsoidal gaussian center, gaussian surround'
            visualRFcenterCharacteristicRadiusDegs = obj.rfComputeStruct.targetVisualRFparamsVector(2:3);

        case 'arbitrary center, gaussian surround'
            visualRFcenterCharacteristicRadiusDegs = modelConstants.visualRFcenterCharacteristicRadiusDegs;

        otherwise
            error('Unknown targetVisualMapScheme: ''%s''.', targetVisualMapScheme);
    end


    retinalProfileColor = [1 0.0 0.2];
    visualProfileColor  = [0 1.0 0.6];
    targetProfileColor  = [1 0.8 0];

    figNo = 169;
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [15 15 1900 1180], 'Color', 'none');

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 5, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);

    % The PSF
    ax = subplot('Position', subplotPosVectors(1,1).v);
    maxPSF = max(obj.theVlambdaWeightedPSFData.vLambdaWeightedData(:));
    plotPSF(ax, obj.theVlambdaWeightedPSFData, maxPSF, maxSpatialSupportDegs, eccDegs, subjectID, false, false, '');


    % The retinal RF center cone map
    ax = subplot('Position', subplotPosVectors(1,2).v);
    if (conesNumInRFcenter== 1)
        plotTitle = sprintf('retinal RF center\n(%d input cone)', conesNumInRFcenter);
    else
        plotTitle = sprintf('retinal RF center\n(%d input cones)', conesNumInRFcenter);
    end
    
    maxRetinalRFprofile = 0.5*max(sum(obj.rfComputeStruct.theRetinalRFcenterConeMap,1));


    % First row: center cone maps with a common max center RF
    maxCenterRetinalRF = 0.1*max(obj.rfComputeStruct.theRetinalRFcenterConeMap(:));

    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theRetinalRFcenterConeMap,  ...
        maxCenterRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, true, false, ...
        'profileColor', retinalProfileColor);

    % The achieved visual RF center map
    ax = subplot('Position', subplotPosVectors(1,3).v);
    achievedVisualRFcenterMap = obj.rfComputeStruct.theFittedVisualRFcenterConeMap;
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        achievedVisualRFcenterMap,  ...
        maxCenterRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, sprintf('visual RF center\n'), true, true, ...
        'profileColor', visualProfileColor);

    
    % The target visual RF center cone map
    ax = subplot('Position', subplotPosVectors(1,4).v);
    RcText = sprintf('Rc:%2.2f arcmin', visualRFcenterCharacteristicRadiusDegs*60);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.targetVisualRFcenterMap, ...
        maxCenterRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, sprintf('target visual RF center (%s)\n(best fit Gaussian line weighting function)', RcText), true, true, ...
        'profileColor', targetProfileColor);
    
   
    % Second row: The surround cone maps with a common max surround RF
    maxSurroundRetinalRF = 0.02*maxCenterRetinalRF;


    % The retinal RF surround cone map
    ax = subplot('Position', subplotPosVectors(2,2).v);
    plotTitle = sprintf('retinal RF surround\n(%s)', strrep(obj.targetVisualRFDoGparams.retinalConePoolingModel, 'arbitrary center cone weights,', ''));
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theRetinalRFsurroundConeMap, ...
        maxSurroundRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, true,false, ...
        'profileColor', retinalProfileColor);


    % The achieved visual RF surround map
    ax = subplot('Position', subplotPosVectors(2,3).v);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.theFittedVisualRFsurroundConeMap, ...
        maxSurroundRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, sprintf('visual RF surround\n'), true, true, ...
        'profileColor', visualProfileColor);

    % The target visual RF surround cone map
    ax = subplot('Position', subplotPosVectors(2,4).v);
    surroundLabel = sprintf('Rs/Rc:%2.2f, integrated S/C:%2.2f', obj.targetVisualRFDoGparams.surroundToCenterRcRatio, obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio);
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        obj.rfComputeStruct.targetVisualRFsurroundMap,  ...
        maxSurroundRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, sprintf('target visual RF surround\n(%s)', surroundLabel), true, true, ...
        'profileColor', targetProfileColor);
    
   

    
    % Third row: composite RFs
    
    retinalCompositeMap = obj.rfComputeStruct.theRetinalRFcenterConeMap-obj.rfComputeStruct.theRetinalRFsurroundConeMap;
    maxCompositeRetinalRF = maxSurroundRetinalRF;

    % The retinal RF composite cone map
    ax = subplot('Position', subplotPosVectors(3,2).v);
    SCintegratedRatio = sum(obj.rfComputeStruct.theRetinalRFsurroundConeMap(:))/sum(obj.rfComputeStruct.theRetinalRFcenterConeMap(:));
    plotTitle = sprintf('retinal composite RF\n(integrated S/C: %2.2f)', SCintegratedRatio);
    
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        retinalCompositeMap , ...
        maxCompositeRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, false, ...
        'profileColor', retinalProfileColor, ...
        'mode', 'composite');

    

    % The achieved visual composite RF map
    ax = subplot('Position', subplotPosVectors(3,3).v);
    SCintegratedRatio = sum(obj.rfComputeStruct.theFittedVisualRFsurroundConeMap(:))/sum(obj.rfComputeStruct.theFittedVisualRFcenterConeMap(:));
    plotTitle = sprintf('visual composite RF\n(integrated S/C: %2.2f)', SCintegratedRatio);
    achievedVisualCompositeMap = obj.rfComputeStruct.theFittedVisualRFcenterConeMap-obj.rfComputeStruct.theFittedVisualRFsurroundConeMap;
    %maxNegative = max(abs(achievedVisualCompositeMap(achievedVisualCompositeMap<0)));
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        achievedVisualCompositeMap, ...
        maxCompositeRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, true, ...
        'profileColor', visualProfileColor, ...
        'mode', 'composite');


    % The target visual composite RF map
    ax = subplot('Position', subplotPosVectors(3,4).v);
    SCintegratedRatio = sum(obj.rfComputeStruct.targetVisualRFsurroundMap(:))/sum(obj.rfComputeStruct.targetVisualRFcenterMap(:));
    plotTitle = sprintf('target visual composite RF\n(integrated S/C: %2.2f)', SCintegratedRatio);
    targetVisualCompositeMap = obj.rfComputeStruct.targetVisualRFcenterMap-obj.rfComputeStruct.targetVisualRFsurroundMap;
    plotRF(ax, spatialSupportXdegs, spatialSupportYdegs, ...
        targetVisualCompositeMap, ...
        maxCompositeRetinalRF, maxRetinalRFprofile, ...
        maxSpatialSupportDegs, plotTitle, false, true, ...
        'profileColor', targetProfileColor, ...
        'secondProfile', sum(achievedVisualCompositeMap,1), ...
        'secondProfileColor', visualProfileColor, ...
        'secondProfileLineType', '-', ...
        'secondProfileIsShaded', true, ...
        'mode', 'composite');



     % Compute the retinal STFs
    retinalCompositeMap = obj.rfComputeStruct.theRetinalRFcenterConeMap-obj.rfComputeStruct.theRetinalRFsurroundConeMap;
    [~,retinalSTF.center] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.theRetinalRFcenterConeMap,1));
    [~,retinalSTF.surround] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.theRetinalRFsurroundConeMap,1));
    [~,retinalSTF.composite] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(retinalCompositeMap,1));


    % Compute the visual STFs
    visualCompositeMap = obj.rfComputeStruct.theFittedVisualRFcenterConeMap - obj.rfComputeStruct.theFittedVisualRFsurroundConeMap;
    [spatialFrequencySupport,visualSTF.center] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.theFittedVisualRFcenterConeMap,1));
    [~,visualSTF.surround] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.theFittedVisualRFsurroundConeMap,1));
    [~,visualSTF.composite] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(visualCompositeMap,1));

    % Compute the target STFs
    targetCompositeMap = obj.rfComputeStruct.targetVisualRFcenterMap - obj.rfComputeStruct.targetVisualRFsurroundMap;
    [~,targetSTF.center] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.targetVisualRFcenterMap,1));
    [~,targetSTF.surround] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.rfComputeStruct.targetVisualRFsurroundMap,1));
    [~,targetSTF.composite] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(targetCompositeMap,1));

    [~,thePSFMTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            obj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1),sum(obj.theVlambdaWeightedPSFData.vLambdaWeightedData,1));

    maxSTF = max(retinalSTF.center);

    theSTFLabels = {'retinal', 'visual', 'target DoG (visual)'};
    theSTFs = [...
        retinalSTF.center; ...
        visualSTF.center; ...
        targetSTF.center];
        
    theSTFColors = [...
        retinalProfileColor; ...
        visualProfileColor; ...
        targetProfileColor ];

    % Plot the center STFs
    ax = subplot('Position', subplotPosVectors(1,5).v);
    plotSTFs(ax, spatialFrequencySupport, theSTFs, thePSFMTF, theSTFLabels, theSTFColors, maxSTF, true, false, sprintf('center STFs\n'));

    % Plot the surround STFs
    ax = subplot('Position', subplotPosVectors(2,5).v);
    theSTFs = [...
        retinalSTF.surround; ...
        visualSTF.surround; ...
        targetSTF.surround];
    plotSTFs(ax, spatialFrequencySupport, theSTFs, thePSFMTF, {}, theSTFColors, maxSTF, true, false, sprintf('surround STFs\n'));

    % Plot the composite STFs
    ax = subplot('Position', subplotPosVectors(3,5).v);
    theSTFs = [...
        retinalSTF.composite; ...
        visualSTF.composite; ...
        targetSTF.composite];
    stfRMSE = sqrt(1/(numel(obj.rfComputeStruct.theSTF.fitted))*sum((obj.rfComputeStruct.theSTF.fitted - obj.rfComputeStruct.theSTF.target).^2));
    plotSTFs(ax, spatialFrequencySupport, theSTFs, thePSFMTF, {}, theSTFColors, maxSTF, false, false, sprintf('composite STFs (rmse:%2.2f)\n',stfRMSE));

    
    % The retinal RF surround model parameter values & ranges
    ax = subplot('Position', subplotPosVectors(3,1).v);
    RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, ...
       obj.rfComputeStruct.retinalConePoolingParams, ...
       obj.targetVisualRFDoGparams.retinalConePoolingModel);

end


function plotSTFs(ax, spatialFrequencySupport, theSTFs, thePSFMTF, theSTFLabels, theSTFcolors, maxSTF, noXTickLabel, noYTickLabel, plotTitle)

    
    hold(ax, 'on');
    for iSTF = 1:size(theSTFs,1)
        p(iSTF) = plot(ax, spatialFrequencySupport, theSTFs(iSTF,:)/maxSTF, 'o-', 'LineWidth', 1.5, 'MarkerSize', 10, ...
            'MarkerFaceColor', theSTFcolors(iSTF,:), 'MarkerEdgeColor', theSTFcolors(iSTF,:)*0.5, 'Color', theSTFcolors(iSTF,:)*0.2); 
    end
    
    plot(ax, spatialFrequencySupport, thePSFMTF, '-', 'LineWidth', 4, 'Color', [0.8 0.8 0.8]);
    p(numel(p)+1) = plot(ax, spatialFrequencySupport, thePSFMTF, 'k--', 'LineWidth', 2);
    
    axis(ax, 'square');
    set(ax, 'XScale', 'log', 'XLim', [0.3 300], 'YLim', [0.0 1.0], 'YScale', 'linear', ...
            'XTick', [0.1 0.3 1 3 10 30 100 300], ...
            'XTickLabel', {'0.1', '0.3', '1', '3', '10', '30', '100', '300'}, ...
            'YTick', 0:0.2:1, 'FontSize', 16);
    grid(ax, 'on')
    
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    set(ax.XAxis,'TickDir','both', 'TickLength',[0.02, 0.01]);
    
    grid(ax, 'on'); box(ax, 'off');

    if (noXTickLabel)
    else
        xlabel('spatial frequency (c/deg)');
    end

    if (noYTickLabel)
        set(ax, 'YTickLabel', {});
    else
        ylabel('');
    end
    set(ax, 'Color', 'none')
    
    if (~isempty(theSTFLabels))
        theSTFLabels = cat(2, theSTFLabels, {'MTF'});
        hLegend = legend(ax, p, theSTFLabels, 'Location', 'SouthWest');
        legend boxoff
    end
    
    title(ax, plotTitle, 'Color', [0.3 0.3 0.3]);
end


function plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID, noXTickLabel, noYTickLabel, plotTitle)
    psfZLevels = 0.05:0.1:0.95;
    
    contourf(ax,thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.vLambdaWeightedData/maxPSF, psfZLevels);
    hold on;
    midRow = (size(thePSFData.vLambdaWeightedData,1)-1)/2+1;
    plot(ax, thePSFData.psfSupportXdegs, -maxSpatialSupportDegs*0.95 + 1.95*thePSFData.vLambdaWeightedData(midRow,:)/maxPSF*maxSpatialSupportDegs, 'k-', 'LineWidth', 1.5);
    axis(ax,'image'); axis 'xy';
    
    ticks = ticksForSpatialSupport(maxSpatialSupportDegs);
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
            'XTick', ticks, 'YTick', ticks, 'CLim', [0 1], 'FontSize', 16);

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
        title(ax, sprintf('PSF @ (%2.2f, %2.2f) degs\nsubjID: %d', eccDegs(1), eccDegs(2), testSubjectID), 'Color', [0.3 0.3 0.3]);
    else
        title(ax, plotTitle, 'Color', [0.3 0.3 0.3]);
    end
    set(ax, 'color', 'none')
    colormap(ax,brewermap(1024, 'greys'));
end

function plotRetinalSurroundModel(ax, rfSupportX, maxSpatialSupportDegs)
    
    ticks = ticksForSpatialSupport(maxSpatialSupportDegs);
    axis(ax, 'square');
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', [-0.3 1.01], ...
            'XTick', ticks, 'YTick', -1:0.2:1, 'FontSize', 16);
    
    grid(ax, 'on'); box(ax, 'on');
    xtickangle(ax, 90);
    %legend({'target', 'X', 'Y', 'residual'});
    title('surround model');
    xlabel(ax,'degrees');

end




function plotRF(ax, rfSupportX, rfSupportY, RF,  maxRF, maxRFprofile, maxSpatialSupportDegs, ...
    titleString, noXTickLabels, noYTickLabels, varargin)

    p = inputParser;
    p.addParameter('profileColor', [1 0 0], @isnumeric);
    p.addParameter('secondProfile', [], @isnumeric);
    p.addParameter('secondProfileColor', [0 0 0], @isnumeric);
    p.addParameter('secondProfileLineType', '--');
    p.addParameter('secondProfileIsShaded', false, @islogical);
    p.addParameter('mode', 'component', @(x)(ischar(x)&&(ismember(x, {'component','composite'}))));
    p.parse(varargin{:});
    profileColor = p.Results.profileColor;
    secondProfile = p.Results.secondProfile;
    secondProfileColor = p.Results.secondProfileColor;
    secondProfileLineType = p.Results.secondProfileLineType;
    secondProfileIsShaded = p.Results.secondProfileIsShaded;
    mode = p.Results.mode;

    RFprofile = sum(RF,1);
    imagesc(ax, rfSupportX, rfSupportY, RF/maxRF);
    hold on;
    
    if (strcmp(mode, 'component'))
        scale1 = 0.95;
        scale2 = 1.95;
    else
        scale1 = 0.75;
        scale2 = 1.75;
    end

    
    if (~isempty(secondProfile))
        if (secondProfileIsShaded)
            shadedAreaPlot(ax,rfSupportX, ...
                -maxSpatialSupportDegs*scale1 + scale2*secondProfile/maxRFprofile*maxSpatialSupportDegs, ...
                -maxSpatialSupportDegs*scale1, secondProfileColor, secondProfileColor*0.5, 1.0, 1.5, '-');
        else
            plot(ax, rfSupportX, -maxSpatialSupportDegs*scale1 + scale2*secondProfile/maxRFprofile*maxSpatialSupportDegs, 'k', ...
                'Color', secondProfileColor, 'LineWidth', 1.5, 'LineStyle', secondProfileLineType);
        end

    end
    plot(ax, rfSupportX, -maxSpatialSupportDegs*scale1 + rfSupportX*0, 'k-');
    plot(ax, rfSupportX, -maxSpatialSupportDegs*scale1 + scale2*RFprofile/maxRFprofile*maxSpatialSupportDegs, 'r-', 'Color', profileColor, 'LineWidth', 1.5);
    

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

    if (strcmp('mode', 'component'))
        text(0,0,sprintf(' vol: %2.3f', sum(RF(:))));
    end
    set(ax, 'Color', 'none');
    
    grid(ax, 'on');
    xtickangle(ax, 90);
    colormap(ax,brewermap(1024, 'greys'));
    title(ax, titleString, 'Color', [0.3 0.3 0.3]);
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
