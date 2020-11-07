function visualizeConeWeights(obj)

    xRange = []; yRange = [];
    rgcsNum = numel(obj.rgcRFspacingsDegs);
    visualizedRGCindices = round(rgcsNum/2) +[-5 0 5];
    
    for iRGC = 1:numel(visualizedRGCindices)
        RGCindex = visualizedRGCindices(iRGC);
        % RGC data
        rgcSpacing = obj.rgcRFspacingsDegs(RGCindex);
        rgcPosition = obj.rgcRFpositionsDegs(RGCindex,:);

        % RF center connections data
        connectivityVector  = full(squeeze(obj.coneWeights.center(:, RGCindex)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        center.coneWeights = squeeze(full(obj.coneWeights.center(coneIndicesConnectedToCenter, RGCindex)));
        center.conePositions = obj.inputConeMosaicMetaData.conePositionsDegs(coneIndicesConnectedToCenter,:);
        center.coneSpacings = obj.inputConeMosaicMetaData.coneSpacingsDegs(coneIndicesConnectedToCenter);

         % RF surround connections data
        connectivityVector  = full(squeeze(obj.coneWeights.surround(:, RGCindex)));
        coneIndicesConnectedToSurround = find(connectivityVector>0);
        surround.coneWeights = squeeze(full(obj.coneWeights.surround(coneIndicesConnectedToSurround, RGCindex)));
        surround.conePositions = obj.inputConeMosaicMetaData.conePositionsDegs(coneIndicesConnectedToSurround,:);
        surround.coneSpacings = obj.inputConeMosaicMetaData.coneSpacingsDegs(coneIndicesConnectedToSurround);

        % Plot
        [xRange, yRange] = plotRGCandConeWeights(1000+iRGC, RGCindex, rgcPosition, rgcSpacing, center, surround, ...
            obj.synthesizedRFparams.visual.centerCharacteristicRadiiDegs(RGCindex), ...
            obj.synthesizedRFparams.visual.centerPeakSensitivities(RGCindex), ...
            obj.synthesizedRFparams.visual.surroundCharacteristicRadiiDegs(RGCindex), ...
            obj.synthesizedRFparams.visual.surroundPeakSensitivities(RGCindex), ...
            obj.synthesizedRFparams.retinal.surroundCharacteristicRadiiDegs(RGCindex), ...
            xRange, yRange);
    end
end

function [xRange, yRange] = plotRGCandConeWeights(figNo, RGCindex, rgcPosition, rgcSpacing, center, surround, ...
    deconvModelCenterVisualCharacteristicRadius, ...
    deconvModelCenterVisualPeakSensitivity, ...
    deconvModelSurroundVisualCharacteristicRadius, ...
    deconvModelSurroundVisualPeakSensitivity, ...
    deconvModelSurroundRetinalCharacteristicRadius, xRange, yRange)

    edgeAlpha = 0.0;
    xOutline = cosd(0:60:360);
    yOutline = sind(0:60:360);
    
    xOutline2 = cosd(0:10:360);
    yOutline2 = sind(0:10:360);
    
    if (isempty(xRange))
        xRange = [...
            min(surround.conePositions(:,1))-mean(surround.coneSpacings(:)) ...
            max(surround.conePositions(:,1))+mean(surround.coneSpacings(:))];
        xRange = mean(xRange) + 0.5*diff(xRange)*[-1 1];

        yRange = [...
            min(surround.conePositions(:,2))-mean(surround.coneSpacings(:)) ...
            max(surround.conePositions(:,2))+mean(surround.coneSpacings(:))];
        yRange = xRange - mean(xRange) + mean(yRange);
    end
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1400 700], 'Name', sprintf('RGC #%d', RGCindex));
    ax = subplot('Position', [0.05 0.05 0.40 0.90]);
    hold(ax, 'on');
    
    % The cones feeding into the surround
    edgeColor = [0 0 1];
    faceColor = [0 0 1];
    for iCone = 1:numel(surround.coneWeights)
        coneOutlineX = xOutline * surround.coneSpacings(iCone) * 0.5 + surround.conePositions(iCone,1);
        coneOutlineY = yOutline * surround.coneSpacings(iCone) * 0.5 + surround.conePositions(iCone,2);
        faceAlpha = min([1 3*surround.coneWeights(iCone)/max(center.coneWeights)]);
        patchContour(ax, coneOutlineX, coneOutlineY, faceColor, edgeColor, faceAlpha, edgeAlpha);
    end
    
    % The cones feeding into the center
    edgeColor = [1 0 0];
    faceColor = [1 0 0];
    for iCone = 1:numel(center.coneWeights)
        coneOutlineX = xOutline * center.coneSpacings(iCone) * 0.5 + center.conePositions(iCone,1);
        coneOutlineY = yOutline * center.coneSpacings(iCone) * 0.5 + center.conePositions(iCone,2);
        faceAlpha = min([1 3*center.coneWeights(iCone)/max(center.coneWeights)]);
        patchContour(ax, coneOutlineX, coneOutlineY, faceColor, edgeColor, faceAlpha, edgeAlpha);
    end
    
    % The corresponding visual characteristic radius
    plot(xOutline2 * deconvModelCenterVisualCharacteristicRadius * 0.5 + rgcPosition(1), ...
         yOutline2 * deconvModelCenterVisualCharacteristicRadius * 0.5 + rgcPosition(2), 'k--', 'LineWidth', 1.5);
    
    % The visual RF profile for the center
    gain = 0.95*(yRange(end)-yRange(1)) / deconvModelCenterVisualPeakSensitivity;
    xx = linspace(xRange(1), xRange(end), 100);
    plot(xx,  yRange(1) + gain*deconvModelCenterVisualPeakSensitivity*exp(-((xx-rgcPosition(1))/deconvModelCenterVisualCharacteristicRadius).^2), ...
        'r-', 'LineWidth', 1.5);
    plot(xx,  yRange(1) + gain*deconvModelSurroundVisualPeakSensitivity*exp(-((xx-rgcPosition(1))/deconvModelSurroundVisualCharacteristicRadius).^2), ...
        'b-', 'LineWidth', 1.5);
    
    % The RGC center outline
    edgeColor = [0 0 0];
    faceColor = [0.8 0.8 0.8];
    faceAlpha = 0.0;
    rgcOutlineX = xOutline * rgcSpacing * 0.5 + rgcPosition(1);
    rgcOutlineY = yOutline * rgcSpacing * 0.5 + rgcPosition(2);
    patchContour(ax, rgcOutlineX, rgcOutlineY, faceColor, edgeColor, faceAlpha, edgeAlpha);
    
    % Finish plot
    axis(ax, 'equal');
    set(ax, 'XLim', xRange, 'YLim', yRange);
    
    ax = subplot('Position', [0.55 0.05 0.40 0.90]);
    hold(ax, 'on');
    % The cones feeding into the surround
    edgeColor = [0 0 1];
    faceColor = [0 0 1];
    for iCone = 1:numel(surround.coneWeights)
        coneOutlineX = xOutline * surround.coneSpacings(iCone) * 0.5 + surround.conePositions(iCone,1);
        coneOutlineY = yOutline * surround.coneSpacings(iCone) * 0.5 + surround.conePositions(iCone,2);
        faceAlpha = surround.coneWeights(iCone)/max(surround.coneWeights);
        patchContour(ax, coneOutlineX, coneOutlineY, faceColor, edgeColor, faceAlpha, edgeAlpha);
    end
    
    plot(xOutline2 * deconvModelSurroundRetinalCharacteristicRadius* 0.5 + rgcPosition(1), ...
         yOutline2 * deconvModelSurroundRetinalCharacteristicRadius* 0.5 + rgcPosition(2), 'g-', 'LineWidth', 1.5);
     
    plot(xOutline2 * deconvModelSurroundVisualCharacteristicRadius * 0.5 + rgcPosition(1), ...
         yOutline2 * deconvModelSurroundVisualCharacteristicRadius * 0.5 + rgcPosition(2), 'k--', 'LineWidth', 1.5);
     
    % The RGC center outline
    edgeColor = [0 0 0];
    faceColor = [0.8 0.8 0.8];
    faceAlpha = 0.0;
    rgcOutlineX = xOutline * rgcSpacing * 0.5 + rgcPosition(1);
    rgcOutlineY = yOutline * rgcSpacing * 0.5 + rgcPosition(2);
    patchContour(ax, rgcOutlineX, rgcOutlineY, faceColor, edgeColor, faceAlpha, edgeAlpha);
    
    % Finish plot
    axis(ax, 'equal');
    set(ax, 'XLim', xRange, 'YLim', yRange);
end


function patchContour(theAxes, xOutline, yOutline, faceColor, edgeColor, faceAlpha, edgeAlpha)
    v = [xOutline(:) yOutline(:)];
    f = 1:numel(xOutline);
    patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, ... 
           'EdgeAlpha', edgeAlpha, 'LineWidth', 1.5);
end

