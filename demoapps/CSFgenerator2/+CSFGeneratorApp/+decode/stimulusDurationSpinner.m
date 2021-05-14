function stimulusDurationSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusDurationSpinner.Value = value*1000;
        case 'sliderToValue'
            app.stimParams.durationSec = value/1000;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end