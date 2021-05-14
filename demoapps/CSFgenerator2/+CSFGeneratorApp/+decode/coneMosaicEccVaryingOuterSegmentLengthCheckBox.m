function coneMosaicEccVaryingOuterSegmentLengthCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingOuterSegmentLengthCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingOuterSegmentLength = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end