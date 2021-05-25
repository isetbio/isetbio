function coneMosaicEccVaryingOuterSegmentLengthCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingOuterSegmentLengthCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingOuterSegmentLength = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.eccVaryingOuterSegmentLength = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end