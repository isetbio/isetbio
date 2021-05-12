function roiRadialEccentricityScalingSwitch(app, direction, value)
    
    if (strcmp(value, 'log'))
        CSFGeneratorApp.initialize.radialEccentricityWithLogScaling(app);
    else
        CSFGeneratorApp.initialize.radialEccentricityWithLinearScaling(app);
    end
    
    switch direction
        case 'valueToSlider'
            app.roiRadialEccentricityScalingSwitch.Value = value;
        case 'sliderToValue'
            app.roiParams.radialEccentricityScaling = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
