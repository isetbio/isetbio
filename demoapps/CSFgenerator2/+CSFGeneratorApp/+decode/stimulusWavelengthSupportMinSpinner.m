function stimulusWavelengthSupportMinSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusWavelengthSupportMinSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.wavelengthSupportMin = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end