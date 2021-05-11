function visualFieldView(app, mode)

    switch (mode)
        case 'initialize'
            initializeVisualFieldView(app);
        case 'update'
            updateVisualFieldViewWithNewData(app);
    end
end


function initializeVisualFieldView(app)
    % Plot the crosshairs
    plot(app.visualFieldView, app.visualFieldParams.maxEcc*[-1 1], [0 0], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    hold(app.visualFieldView, 'on');
    plot(app.visualFieldView, [0 0], app.visualFieldParams.maxEcc*[-1 1], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    % Plot the eccentricity rings
    eccRadii = app.visualFieldRadialEccentricitySlider.MajorTickLabels;
    for ir = 1:numel(eccRadii)
        logRadius = ir/numel(eccRadii)*app.visualFieldParams.maxEcc;
        plot(app.visualFieldView, logRadius*cosd(0:5:360), logRadius*sind(0:5:360), 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
            
    % Plot the optic disk location
    odOutline.x = nan;
    odOutline.y = nan;
    odColor = [1 0.3 0];
    app.odOnVisualFieldPlotHandle = patch(app.visualFieldView, odOutline.x, odOutline.y, odColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    % Plot the mosaic outline
    mosaicOutline.x = nan;
    mosaicOutline.y = nan;
    app.mosaicOnVisualFieldPlotHandle = patch(app.visualFieldView, mosaicOutline.x, mosaicOutline.y, [0.1 0.95 0.7], 'FaceAlpha', 0.5, 'EdgeColor', [0.1 0.3 0.8]);   % scatter(app.visualFieldView, 0, 0, 14*14, 'ro', 'MarkerEdgeColor', [0.1 0.3 0.8], 'MarkerFaceColor', [0.1 0.95 0.7], 'MarkerFaceAlpha', 0.5);
            
    hold(app.visualFieldView, 'off');
    set(app.visualFieldView, 'XTick', 0, 'YTick', 0, 'XLim', app.visualFieldParams.maxEcc*[-1 1], 'YLim', app.visualFieldParams.maxEcc*[-1 1], 'FontSize', 14);
    set(app.visualFieldView, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    set(app.visualFieldView, 'XTickLabel', '0');
    set(app.visualFieldView, 'YTickLabel', '0');
    box(app.visualFieldView, 'off');
    xlabel(app.visualFieldView, 'x (degs)');
    ylabel(app.visualFieldView, 'y (degs)');
            
    % Do not show the interactions toolbax
    app.visualFieldView.Toolbar.Visible = 'off';
            
    % No interactions
    app.visualFieldView.Interactions = [];     
end

function updateVisualFieldViewWithNewData(app)

    % Update the OD
    [xx,yy] = generateOpticDiskOutline(app);
    set(app.odOnVisualFieldPlotHandle, 'XData', xx, 'YData', yy);
        
    % Update the mosaic outline
    mosaicOutline = generateMosaicOutlineForRetinalPatchView(app);
    set(app.mosaicOnVisualFieldPlotHandle, ...
                    'XData', mosaicOutline.x, ...
                    'YData', mosaicOutline.y);
                
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
    r = sqrt(x.^2+y.^2);
    r = app.visualFieldParams.maxEcc * log10(1+r)/log10(1+app.visualFieldParams.maxEcc);
    theta = atan2(y,x);
    xx = r .* cos(theta);
    yy = r .* sin(theta);
end

function mosaicOutline = generateMosaicOutlineForRetinalPatchView(app)
    x = app.coneMosaicParams.eccentricityDegs(1) + 0.5*app.coneMosaicParams.sizeDegs(1) * [-1 -1 1 1 -1];
    y = app.coneMosaicParams.eccentricityDegs(2) + 0.5*app.coneMosaicParams.sizeDegs(2) * [-1 1 1 -1 -1];

    r = sqrt(x.^2+y.^2);
    r = app.visualFieldParams.maxEcc * log10(1+r)/log10(1+app.visualFieldParams.maxEcc);
    theta = atan2(y,x);
    mosaicOutline.x = r .* cos(theta);
    mosaicOutline.y = r .* sin(theta);
end
        