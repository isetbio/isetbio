function roiPolarEccentricitySlider(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.roiPolarEccentricityKnob.Value = 360-value;
        case 'sliderToValue'
            app.roiParams.polarEccentricityDegs = 360-value;
             
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update the text
    app.polareccentricitydegsLabel.Text = sprintf('polar eccentricity (degs): %2.0f', app.roiParams.polarEccentricityDegs);
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.eccentricityDegs(1) = ...
        app.roiParams.radialEccentricityDegs * cosd(app.roiParams.polarEccentricityDegs);
    app.coneMosaicParams.eccentricityDegs(2) = ...
        app.roiParams.radialEccentricityDegs * sind(app.roiParams.polarEccentricityDegs);       
    
    % Update the x,y eccentricity on the GUI
    app.roiEccentricityX.Value = app.coneMosaicParams.eccentricityDegs(1);
    app.roiEccentricityY.Value = app.coneMosaicParams.eccentricityDegs(2);
end

