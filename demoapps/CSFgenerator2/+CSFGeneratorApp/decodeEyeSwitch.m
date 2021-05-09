function decodeEyeSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.visualFieldEyeSwitch.Value = value;
        case 'sliderToValue'
            app.visualFieldParams.whichEye = value;
            
        otherwise
            error('Unknown FieldOfViewSlider.direction; ''%s''.'\n', direction);
    end
    
    % Set the associated coneMosaicParam
    app.coneMosaicParams.whichEye = value;
end
