function roiEyeSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.roiEyeSwitch.Value = value;
        case 'sliderToValue'
            app.roiParams.whichEye = value;
            
        otherwise
            error('Unknown FieldOfViewSlider.direction; ''%s''.'\n', direction);
    end
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.whichEye = value;
end
