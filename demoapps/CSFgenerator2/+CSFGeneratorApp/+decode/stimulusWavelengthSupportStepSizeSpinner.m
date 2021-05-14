function stimulusWavelengthSupportStepSizeSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusWavelengthSupportStepSizeSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.wavelengthSupportStepSize = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end