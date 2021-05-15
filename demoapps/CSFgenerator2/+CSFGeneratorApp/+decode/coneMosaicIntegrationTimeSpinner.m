function coneMosaicIntegrationTimeSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicIntegrationTimeSpinner.Value = value*1000;
        case 'sliderToValue'
            app.coneMosaicParams.integrationTimeSeconds = value/1000;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.integrationTime = value/1000;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end