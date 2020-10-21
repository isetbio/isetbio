function coneDensity2DMaps()

    % Spatial sampling of the retina in microns
    conePosXmicrons = -6000:100:6000;
    conePosYmicrons = -6000:100:6000;
    [X,Y] = meshgrid(conePosXmicrons, conePosYmicrons);

    conePosMicrons = zeros(numel(X),2);
    conePosMicrons(:,1) = X(:);
    conePosMicrons(:,2) = Y(:);
    
    % Left eye
    whichEye = RGCmodels.Watson.constants.leftEye;
    [~, ~, leftRetinaConeDensityMapDegs2] = RGCmodels.Watson.compute.coneSpacingAtRetinalPositions(whichEye, conePosMicrons);
    leftRetinaConeDensityMapDegs2 = reshape(leftRetinaConeDensityMapDegs2, size(X));
    
    % Right eye
    whichEye = RGCmodels.Watson.constants.rightEye;
    [~, ~, rightRetinaConeDensityMapDegs2] = RGCmodels.Watson.compute.coneSpacingAtRetinalPositions(whichEye, conePosMicrons);
    rightRetinaConeDensityMapDegs2 = reshape(rightRetinaConeDensityMapDegs2, size(X));
    
    % Plot
    figure(2);
    ax = subplot(1,2,1);
    plotData(ax, conePosXmicrons, conePosYmicrons, leftRetinaConeDensityMapDegs2, [100 20000], 'left retina');
    
    ax = subplot(1,2,2);
    plotData(ax, conePosXmicrons, conePosYmicrons, rightRetinaConeDensityMapDegs2, [100 20000], 'right retina');
end

function plotData(ax, conePosXmicrons, conePosYmicrons, densityMap, densityRange, titleLabel)
    contourf(ax,conePosXmicrons, conePosYmicrons, log10(densityMap), linspace(log10(densityRange(1)), log10(densityRange(2)),20));
    hold(ax, 'on');
    m1 = min(densityMap(:));
    [m2, idx] = max(densityMap(:));
    [peakRow, peakCol] = ind2sub(size(densityMap), idx);
    xMin = min(conePosXmicrons);
    xMax = max(conePosXmicrons);
    yMin = min(conePosYmicrons);
    yMax = max(conePosYmicrons);
    plot(ax,conePosXmicrons, yMin+10*(yMax-yMin)*(densityMap(peakRow,:)-m1)/(m2-m1), 'k-', 'LineWidth', 1.5);
    axis(ax,'equal')
    set(ax, 'FontSize', 16, 'XLim', [xMin, xMax], 'YLim', [yMin, yMax], 'CLim', [log10(densityRange(1)) log10(densityRange(2))]);
    xlabel(ax, 'retinal microns');
    title(ax,titleLabel);
end

