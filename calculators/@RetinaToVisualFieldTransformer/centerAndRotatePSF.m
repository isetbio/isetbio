function data = centerAndRotatePSF(data)
    % Generate a binary image
    binaryImage = imbinarize(data/max(data(:)));

    % Determine centroid and orientation of major axis with x-axis
    measurements = regionprops(binaryImage, data, {'WeightedCentroid', 'Orientation'});
    weightedCentroid = measurements.WeightedCentroid;
    orientation = measurements.Orientation;

    % Determine translation vector to center the data
    center = [size(data,2) size(data,1)]/2;
    translationVector = center - weightedCentroid;

    % Center the data
    data = imtranslate(data, translationVector, 'FillValues', 0);

    % Unrotate the data
    data = imrotate(data, -orientation, 'bilinear', 'crop');
end