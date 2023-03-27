function mosaicParams = getMosaicParams(mosaicEcc)
    % Method to obtain the (x,y) eccentricity and (x,y) size of an mRGC
    % mosaic based on some descriptor (for now 1D mosaicEcc)

    switch (mosaicEcc)
        case 0
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
            % which covers the [1 - 4] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [0 0], ...
                'sizeDegs', [3 3]);

        case 2.5
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
            % which covers the [1 - 4] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [2.5 0], ...
                'sizeDegs', [3 3]);

        case 7.0
            % Mosaic params to employ. This is for the 7.0 deg - centered mosaic
            % which covers the [4-10] deg eccentricity range
            mosaicParams = struct(...
                'eccDegs', [7 0], ...
                'sizeDegs', [6 3]);

        case -10.0
            mosaicParams = struct(...
                'eccDegs', [-10 0], ...
                'sizeDegs', [6 3]);

        otherwise
            error('No data for this eccentricity')
    end
end
