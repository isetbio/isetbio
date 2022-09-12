function theCircularPSF = circularlySymmetricPSF(thePSF, mode)
    
    switch (mode)
        case RetinaToVisualFieldTransformer.psfCircularSymmetryModeNone
            theCircularPSF = thePSF;
            
        case RetinaToVisualFieldTransformer.psfCircularSymmetryModeAverage
            theCircularPSF = circularAverage(thePSF);

        case RetinaToVisualFieldTransformer.psfCircularSymmetryModeBestResolution
            thePSF =  centerAndRotatePSF(thePSF);
            theCircularPSF = circularFromSlice(thePSF,2);

        case RetinaToVisualFieldTransformer.psfCircularSymmetryModeBestResolution
            thePSF =  centerAndRotatePSF(thePSF);
            theCircularPSF = circularFromSlice(thePSF,1);

        otherwise
            error('Unknown circular symmetric PSF mode: ''%s''.', mode);
    end
end

function dataOut = circularFromSlice(dataIn, sliceDimension)

    switch (sliceDimension)
        case 1
            midY = (size(dataIn,1)-1)/2+1;
            slice = squeeze(dataIn(midY,:));
            slice = slice(:);
        case 2
            midX = (size(dataIn,2)-1)/2+1;
            slice = squeeze(dataIn(:,midX));
            slice = slice(:);
    end

    
    n = numel(slice);
    x = 1:n; x = x - mean(x);

    [X,Y] = meshgrid(x,x);
    Rdest = sqrt(X.^2+Y.^2);

    Rsource = abs(x);

    dataOut = zeros(n,n);

    radii = 0:1:max(x);
    for q = 1:(length(radii)-1)
        indexDest = find( ...
            (Rdest >= radii(q)) & ...
            (Rdest <  radii(q + 1)) ...
            );
        indexSource = find( ...
            (Rsource >= radii(q)) & ...
            (Rsource <  radii(q + 1)) ...
            );
        weights = 1.0 - (Rsource(indexSource)-radii(q));
        if (~isempty(indexSource)) && (~isempty(indexDest))
            dataOut(indexDest) = sum(weights .* (reshape( slice(indexSource) , size(weights)) ),2) / sum(weights,2); 
        end
    end

    dataOut = dataOut / sum(dataOut(:)) * sum(dataIn(:));

end


function dataOut = circularAverage(dataIn)
    [m, n] = size(dataIn);
    if (n ~= m)
        error('Input must be a square matrix');
    end

    x = 1:n; x = x - mean(x);

    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2+Y.^2);

    dataOut = zeros(m,n);

    radii = 0:1:max(x);
    for q = 1:(length(radii)-1)
        index = find( ...
            (R >= radii(q)) & ...
            (R <  radii(q + 1)) ...
            );
        weights = 1.0 - (R(index)-radii(q));
        if (~isempty(index))
            dataOut(index) = sum(weights .* dataIn(index), 1) / sum(weights(:)); 
        end
    end

    dataOut = dataOut / sum(dataOut(:)) * sum(dataIn(:));

end

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
