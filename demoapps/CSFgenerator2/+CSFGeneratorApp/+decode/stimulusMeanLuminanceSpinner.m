function stimulusMeanLuminanceSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusMeanLuminanceSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.meanLuminanceCdM2 = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end