function coneMosaicEccVaryingConeApertureBlurCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingConeApertureBlurCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingConeApertureBlur = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.eccVaryingConeBlur = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end