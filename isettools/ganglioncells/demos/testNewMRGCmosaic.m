function testNewMRGCmosaic()

    % Eccentricity
    eccDegs = [12.5 -3.4];  % near the optic disk
    sizeDegs = 2*[1 1];

    % Default creation
    m = mRGCMosaic('eccentricityDegs', eccDegs, 'sizeDegs', sizeDegs);
    m.whichEye
    computeStats(m);
    visualizeMosaicAndItsInput(1,m);

    % Create given an input cone mosaic
    c = cMosaic('eccentricityDegs', eccDegs, 'sizeDegs', m.inputConeMosaic.sizeDegs);
    m2 = mRGCMosaic('inputConeMosaic',c);
    m2.whichEye
    
    computeStats(m);
    visualizeMosaicAndItsInput(2,m2);
end

function computeStats(m)
% Compute densities within confines of mRGC mosaic
    minEcc = m.eccentricityDegs - m.sizeDegs/2;
    maxEcc = m.eccentricityDegs + m.sizeDegs/2;
    rgcsNum = find( ...
        m.rgcRFpositionsDegs(:,1) >= minEcc(1) & ...
        m.rgcRFpositionsDegs(:,1) <= maxEcc(1) & ...
        m.rgcRFpositionsDegs(:,2) >= minEcc(2) & ...
        m.rgcRFpositionsDegs(:,2) <= maxEcc(2));
    LconesNum = find(...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.lConeIndices,1) >= minEcc(1) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.lConeIndices,1) <= maxEcc(1) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.lConeIndices,2) >= minEcc(2) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.lConeIndices,2) <= maxEcc(2));
    MconesNum = find(...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.mConeIndices,1) >= minEcc(1) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.mConeIndices,1) <= maxEcc(1) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.mConeIndices,2) >= minEcc(2) & ...
        m.inputConeMosaic.coneRFpositionsDegs(m.inputConeMosaic.mConeIndices,2) <= maxEcc(2));
    areaDegs2 = (maxEcc(1)-minEcc(1)) * (maxEcc(2)-minEcc(2));
    rgcMeanDensity = numel(rgcsNum)/areaDegs2
    coneMeanDensity = (numel(LconesNum)+numel(MconesNum))/areaDegs2
    mRGCsToConesRatio = rgcMeanDensity/coneMeanDensity
end

function visualizeMosaicAndItsInput(figNo,m)
    dTheta = 30;
    xOutline = 0.5*cosd(0:dTheta:360);
    yOutline = 0.5*sind(0:dTheta:360);

    figure(figNo);clf;
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    hold(ax, 'on');
    

    coneIndices = m.inputConeMosaic.lConeIndices;
    for k = 1:numel(coneIndices)
        kk = coneIndices(k);
        plotRFoutline(ax, ...
            m.inputConeMosaic.coneRFpositionsDegs(kk,1), m.inputConeMosaic.coneRFpositionsDegs(kk,2), ...
            m.inputConeMosaic.coneRFspacingsDegs(kk), ...
            xOutline, yOutline, [1 0 0], 1.0);
    end

    coneIndices = m.inputConeMosaic.mConeIndices;
    for k = 1:numel(coneIndices)
        kk = coneIndices(k);
        plotRFoutline(ax, ...
            m.inputConeMosaic.coneRFpositionsDegs(kk,1), m.inputConeMosaic.coneRFpositionsDegs(kk,2), ...
            m.inputConeMosaic.coneRFspacingsDegs(kk), ...
            xOutline, yOutline, [0 1 0], 1.0);
    end

    coneIndices = m.inputConeMosaic.sConeIndices;
    for k = 1:numel(coneIndices)
        kk = coneIndices(k);
        plotRFoutline(ax, ...
            m.inputConeMosaic.coneRFpositionsDegs(kk,1), m.inputConeMosaic.coneRFpositionsDegs(kk,2), ...
            m.inputConeMosaic.coneRFspacingsDegs(kk), ...
            xOutline, yOutline, [0 0 1], 1.0);
    end

    for k = 1:size(m.rgcRFpositionsDegs,1)
        plotRFoutline(ax, ...
            m.rgcRFpositionsDegs(k,1), m.rgcRFpositionsDegs(k,2), m.rgcRFspacingsDegs(k), ...
            xOutline, yOutline, [0.6 0.6 0.3], 0.0);
    end

    axis(ax,'equal')
end

function plotRFoutline(ax,xo,yo,radius,xOutline,yOutline,color,lineWidth)
    xx = xo + radius*xOutline;
    yy = yo + radius*yOutline;
    if (lineWidth > 0)
        plot(ax, xx, yy, '-', 'Color', color, 'LineWidth', lineWidth);
    else
        patch(xx, yy, [0 0 0], 'EdgeColor', color*0.5, 'FaceAlpha', 0.5, ...
        'FaceColor', color, 'LineWidth', 1.0, ...
        'LineStyle', '-', ...
        'Parent', ax);
    end
end
