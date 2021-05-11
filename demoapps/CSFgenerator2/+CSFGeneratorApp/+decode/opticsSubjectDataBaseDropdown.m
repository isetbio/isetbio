function opticsSubjectDataBaseDropdown(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.opticsSubjectDataBaseDropDown.Value = value;
        case 'sliderToValue'
            app.opticsParams.subjectDataset = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
