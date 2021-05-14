function roiView(app, mode)

    switch (mode)
        case 'initialize'
            initializeROIView(app);
        case 'update'
            updateROIViewWithNewData(app);
    end
end

function initializeROIView(app)
    % Plot the crosshairs
    plot(app.roiView, app.roiParams.maxEcc*2*[-1 1], [0 0], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    hold(app.roiView, 'on');
    plot(app.roiView, [0 0], app.roiParams.maxEcc*2*[-1 1], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    % Plot the eccentricity rings
    eccRadii = app.roiRadialEccentricitySlider.MajorTickLabels;
    for ir = 1:numel(eccRadii)
        radius = ir/numel(eccRadii)*app.roiParams.maxEcc;
        plot(app.roiView, radius*cosd(0:5:360), radius*sind(0:5:360), ...
            'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
            
    % Plot the optic disk location
    odOutline.x = nan;
    odOutline.y = nan;
    odColor = [1 0.5 0.2];
    app.opticDiskOnROIPlotHandle = patch(app.roiView, odOutline.x, odOutline.y, odColor, ...
        'FaceAlpha', 0.5, 'EdgeColor', odColor*0.8, 'LineWidth', 1.0);

    % Plot the mosaic outline
    mosaicOutline.x = nan;
    mosaicOutline.y = nan;
    mosaicColor = [0.1 0.95 0.7];
    app.mosaicOnROIPlotHandle = patch(app.roiView, mosaicOutline.x, mosaicOutline.y, ...
        mosaicColor, 'FaceAlpha', 0.5, 'EdgeColor', mosaicColor*0.8, 'LineWidth', 1.0);
            
    % Plot the stimulus outline
    stimulusOutline.x = nan;
    stimulusOutline.y = nan;
    stimulusColor = [0.8 0.2 0.9];
    app.stimulusOnROIPlotHandle = patch(app.roiView, stimulusOutline.x, stimulusOutline.y, ...
        stimulusColor, 'FaceAlpha', 0.5, 'EdgeColor', stimulusColor*0.8, 'LineWidth', 1.0);
    
    text(app.roiView, app.roiParams.maxEcc-7, app.roiParams.maxEcc, 'cone mosaic', ...
        'FontSize', 16, 'FontSize', 16, 'Color', mosaicColor*0.8);
    text(app.roiView, app.roiParams.maxEcc-4, app.roiParams.maxEcc-3, 'stimulus', ...
        'FontSize', 16, 'FontSize', 16, 'Color', stimulusColor*0.8);
    text(app.roiView, app.roiParams.maxEcc-5, app.roiParams.maxEcc-6, 'optic disk', ...
        'FontSize', 16, 'FontSize', 16, 'Color',odColor*0.8);
    
    hold(app.roiView, 'off');
    set(app.roiView, 'XTick', 0, 'YTick', 0, ...
        'XLim', (app.roiParams.maxEcc+2)*[-1 1], 'YLim', (app.roiParams.maxEcc+2)*[-1 1], 'FontSize', 1);
    set(app.roiView, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    set(app.roiView, 'XTickLabel', '0');
    set(app.roiView, 'YTickLabel', '0');
    box(app.roiView, 'off');
    xlabel(app.roiView, 'x (degs)');
    ylabel(app.roiView, 'y (degs)');
            
    % Do not show the interactions toolbax
    app.roiView.Toolbar.Visible = 'off';
            
    % No interactions
    app.roiView.Interactions = [];     
end

function updateROIViewWithNewData(app)

    % Update the OD
    [xx,yy] = generateOpticDiskOutline(app);
    set(app.opticDiskOnROIPlotHandle, 'XData', xx, 'YData', yy);
        
    % Update the mosaic outline
    mosaicOutline = generateMosaicOutlineForRetinalPatchView(app);
    set(app.mosaicOnROIPlotHandle, ...
                    'XData', mosaicOutline.x, ...
                    'YData', mosaicOutline.y);
                
    % Update the stimulus outline
    stimulusOutline = generateStimulusOutlineForRetinalPatchView(app);
    set(app.stimulusOnROIPlotHandle, ...
                    'XData', stimulusOutline.x, ...
                    'YData', stimulusOutline.y);
end


function [xx,yy] = generateOpticDiskOutline(app)
    [~,odStructDegs] = app.components.coneMosaic.odStruct();
    a = 0.5*odStructDegs.minorAxisDiameter; 
    b = 0.5*odStructDegs.majorAxisDiameter; 
    c = odStructDegs.center;
    cosOutline = cosd(0:10:360);
    sinOutline = sind(0:10:360);
    x = c(1) + a*cosOutline*cosd(odStructDegs.rotation) - b*sinOutline*sind(odStructDegs.rotation);
    y = c(2) + b*sinOutline*cosd(odStructDegs.rotation) + a*cosOutline*sind(odStructDegs.rotation);
    if (strcmp(app.roiParams.radialEccentricityScaling, 'log'))
        r = sqrt(x.^2+y.^2);
        r = app.roiParams.maxEcc * log10(1+r)/log10(1+app.roiParams.maxEcc);
        theta = atan2(y,x);
        xx = r .* cos(theta);
        yy = r .* sin(theta);
    else
        xx = x;
        yy = y;
    end
    
end

function mosaicOutline = generateMosaicOutlineForRetinalPatchView(app)
    x = app.coneMosaicParams.eccentricityDegs(1) + 0.5*app.coneMosaicParams.sizeDegs(1) * [-1 -1 1 1 -1];
    y = app.coneMosaicParams.eccentricityDegs(2) + 0.5*app.coneMosaicParams.sizeDegs(2) * [-1 1 1 -1 -1];

    if (strcmp(app.roiParams.radialEccentricityScaling, 'log'))
        r = sqrt(x.^2+y.^2);
        r = app.roiParams.maxEcc * log10(1+r)/log10(1+app.roiParams.maxEcc);
        theta = atan2(y,x);
        mosaicOutline.x = r .* cos(theta);
        mosaicOutline.y = r .* sin(theta);
    else
        mosaicOutline.x = x;
        mosaicOutline.y = y;
    end
end
        
function stimulusOutline = generateStimulusOutlineForRetinalPatchView(app)

    if (app.stimParams.mosaicCenteredPosition)
        xo = app.coneMosaicParams.eccentricityDegs(1);
        yo = app.coneMosaicParams.eccentricityDegs(2);
    else
        xo = app.stimParams.positionDegs(1);
        yo = app.stimParams.positionDegs(2);
    end
    
    xx = [];
    yy = [];
    pointsPerSide = 5;
    xx = cat(2, xx, repmat(-1, [1 pointsPerSide]));
    xx = cat(2, xx, repmat(-1, [1 pointsPerSide]));
    xx = cat(2, xx, repmat( 1, [1 pointsPerSide]));
    xx = cat(2, xx, repmat( 1, [1 pointsPerSide]));
    xx = cat(2, xx, repmat(-1, [1 pointsPerSide]));
    
    yy = cat(2, yy, repmat(-1, [1 pointsPerSide]));
    yy = cat(2, yy, repmat( 1, [1 pointsPerSide]));
    yy = cat(2, yy, repmat( 1, [1 pointsPerSide]));
    yy = cat(2, yy, repmat(-1, [1 pointsPerSide]));
    yy = cat(2, yy, repmat(-1, [1 pointsPerSide]));
    
    x = xo + 0.5*app.stimParams.sizeDegs * xx;
    y = yo + 0.5*app.stimParams.sizeDegs * yy;

    if (strcmp(app.roiParams.radialEccentricityScaling, 'log'))
        r = sqrt(x.^2+y.^2);
        r = app.roiParams.maxEcc * log10(1+r)/log10(1+app.roiParams.maxEcc);
        theta = atan2(y,x);
        stimulusOutline.x = r .* cos(theta);
        stimulusOutline.y = r .* sin(theta);
    else
        stimulusOutline.x = x;
        stimulusOutline.y = y;
    end
end
