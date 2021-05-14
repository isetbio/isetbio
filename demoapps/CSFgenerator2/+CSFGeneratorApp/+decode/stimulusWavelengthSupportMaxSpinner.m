function stimulusWavelengthSupportMaxSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusWavelengthSupportMaxSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.wavelengthSupportMax = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end