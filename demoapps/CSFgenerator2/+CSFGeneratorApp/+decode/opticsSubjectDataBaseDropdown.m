function opticsSubjectDataBaseDropdown(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.opticsSubjectDataBaseDropDown.Value = value;
        case 'sliderToValue'
            app.opticsParams.subjectDataset = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    if (strcmp(value, 'Polans2015'))
        % Only right eye data exist, so disable eye switch
        app.roiEyeSwitch.Value = 'right eye';
        app.roiEyeSwitch.Enable = 'off';
    else
        % Enable eye switch
        app.roiEyeSwitch.Enable = 'on';
    end
end
