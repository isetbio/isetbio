function polarEccentricitySlider(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.visualFieldPolarEccentricityKnob.Value = 360-value;
        case 'sliderToValue'
            app.visualFieldParams.polarEccentricityDegs = 360-value;
             
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update the text
    app.polareccentricitydegsLabel.Text = sprintf('polar eccentricity (degs): %2.0f', app.visualFieldParams.polarEccentricityDegs);
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.eccentricityDegs(1) = ...
        app.visualFieldParams.radialEccentricityDegs * cosd(app.visualFieldParams.polarEccentricityDegs);
    app.coneMosaicParams.eccentricityDegs(2) = ...
        app.visualFieldParams.radialEccentricityDegs * sind(app.visualFieldParams.polarEccentricityDegs);       
end

