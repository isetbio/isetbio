function radialEccentricitySlider(app, direction, value)

    tickLabels = app.visualFieldRadialEccentricitySlider.MajorTickLabels;
    actualValues = zeros(1, numel(tickLabels));
    for k = 1:numel(tickLabels)
        actualValues(k) = str2double(tickLabels{k});
    end
    
    switch direction
        case 'valueToSlider'
            app.visualFieldRadialEccentricitySlider.Value = ...
                interp1(actualValues, 0:(numel(actualValues)-1),value);
        case 'sliderToValue'
            app.visualFieldParams.radialEccentricityDegs = ...
                interp1(0:(numel(actualValues)-1), actualValues, value);
        otherwise
            error('Unknown decodeRadialEccentricitySlider.direction; ''%s''.'\n', direction);
    end
    
    % Update the text
    if (app.visualFieldParams.radialEccentricityDegs < 1.0)
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.2f', app.visualFieldParams.radialEccentricityDegs);
    elseif (app.visualFieldParams.radialEccentricityDegs < 10)
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.1f', app.visualFieldParams.radialEccentricityDegs);
    else
        app.radialeccentricitydegsLabel.Text = sprintf('radial eccentricity (degs): %2.0f', app.visualFieldParams.radialEccentricityDegs);
    end
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.eccentricityDegs(1) = ...
        app.visualFieldParams.radialEccentricityDegs * cosd(app.visualFieldParams.polarEccentricityDegs);
    app.coneMosaicParams.eccentricityDegs(2) = ...
        app.visualFieldParams.radialEccentricityDegs * sind(app.visualFieldParams.polarEccentricityDegs);       
end