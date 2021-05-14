function opticsSubjectIDSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsSubjectIDSpinner.Value = value;
        case 'sliderToValue'
            app.opticsParams.subjectID = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
