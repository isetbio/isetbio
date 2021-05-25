function roiRadialEccentricitySlider(app, direction, value)

    tickLabels = app.roiRadialEccentricitySlider.MajorTickLabels;
    actualValues = zeros(1, numel(tickLabels));
    for k = 1:numel(tickLabels)
        actualValues(k) = str2double(tickLabels{k});
    end
    
    switch direction
        case 'valueToSlider'
            app.roiRadialEccentricitySlider.Value = ...
                interp1(actualValues, 0:(numel(actualValues)-1),value);
        case 'sliderToValue'
            app.roiParams.radialEccentricityDegs = ...
                interp1(0:(numel(actualValues)-1), actualValues, value);
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update the text
    if (app.roiParams.radialEccentricityDegs < 1.0)
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.2f', app.roiParams.radialEccentricityDegs);
    elseif (app.roiParams.radialEccentricityDegs < 10)
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.1f', app.roiParams.radialEccentricityDegs);
    else
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.0f', app.roiParams.radialEccentricityDegs);
    end
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.eccentricityDegs(1) = ...
        app.roiParams.radialEccentricityDegs * cosd(app.roiParams.polarEccentricityDegs);
    app.coneMosaicParams.eccentricityDegs(2) = ...
        app.roiParams.radialEccentricityDegs * sind(app.roiParams.polarEccentricityDegs);  
    
    % Update the x,y eccentricity on the GUI
    app.roiEccentricityX.Value = app.coneMosaicParams.eccentricityDegs(1);
    app.roiEccentricityY.Value = app.coneMosaicParams.eccentricityDegs(2);
end