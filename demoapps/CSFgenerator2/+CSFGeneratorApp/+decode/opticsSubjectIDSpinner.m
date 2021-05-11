function opticsSubjectIDSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsSubjectIDSpinner.Value = value;
        case 'sliderToValue'
            app.opticsParams.subjectID = value;
            % Compute new optics
            CSFGeneratorApp.generate.optics(app);
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
