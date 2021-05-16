function csfSourceSignalKnob(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.csfSourceSignalKnob.Value = value;
        case 'sliderToValue'
            app.csfParams.sourceSignal = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
end
