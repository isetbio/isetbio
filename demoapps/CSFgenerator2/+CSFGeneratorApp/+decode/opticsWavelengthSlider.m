function opticsWavelengthSlider(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsVisualizedWavelengthSlider.Value = value;
        case 'sliderToValue'
            app.opticsParams.visualizedWavelength = value;
            
        otherwise
            error('Unknown opticsWavelengthSlider.direction; ''%s''.'\n', direction);
    end
end

