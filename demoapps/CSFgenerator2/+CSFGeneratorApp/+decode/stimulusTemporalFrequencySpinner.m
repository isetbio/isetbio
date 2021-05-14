function stimulusTemporalFrequencySpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusTemporalFrequencySpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.temporalFrequencyHz = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end