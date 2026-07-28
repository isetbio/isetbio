function theFittedGaussianEllipsoids = analyzeSpatialRFmaps(theRFmaps, spatialSupportDegs, GaussianEllipsoidSamplesNum)

	theFittedGaussianEllipsoids = cell(1, numel(theRFmaps));

    parfor iRGC = 1:numel(theRFmaps)
    	fprintf('Fitting Gaussian ellispoid to RF map %d of %d\n', iRGC, numel(theRFmaps));
        theRFmap = theRFmaps{iRGC};
        rfMax = max(abs(theRFmap(:)));

        assert(...
        	(size(theRFmap,1) == numel(spatialSupportDegs)) && (size(theRFmap,2) == numel(spatialSupportDegs)), ...
        	sprintf('size error! RFmap is %d x %d, spatialSupport is %d', size(theRFmap,1), size(theRFmap,2), numel(spatialSupportDegs)));
        
        % Estimate the RF centroid
        theCentroid = rfCentroid(theRFmap, spatialSupportDegs);

        [theFittedGaussianEllipsoidParams, ...
         theFittedGaussianEllipsoid, ...
         theFittedGaussianEllipsoidSpatialSupportXDegs, ...
         theFittedGaussianEllipsoidSpatialSupportYDegs] = RGCMosaicConstructor.helper.fit.GaussianEllipsoidTo2DRFmap(...
                spatialSupportDegs, spatialSupportDegs, theRFmap/rfMax, GaussianEllipsoidSamplesNum, ...
                'fixedCentroid', theCentroid);

        minorSigmaDegs = min([theFittedGaussianEllipsoidParams.xSigma theFittedGaussianEllipsoidParams.ySigma]);
        majorSigmaDegs = max([theFittedGaussianEllipsoidParams.xSigma theFittedGaussianEllipsoidParams.ySigma]);

        theFittedGaussianEllipsoids{iRGC} = struct(...
            'rfMap', theFittedGaussianEllipsoid * rfMax, ...
            'xSupportDegs', theFittedGaussianEllipsoidSpatialSupportXDegs, ...
            'ySupportDegs', theFittedGaussianEllipsoidSpatialSupportYDegs, ...
            'minorSigmaDegs', minorSigmaDegs, ...
            'majorSigmaDegs', majorSigmaDegs ...
        );
    end
end

function theCentroid = rfCentroid(theRFmap, spatialSupportDegs)
	blurFilter = fspecial('gaussian', 10, 2);
	theBlurredRFmap = conv2(theRFmap, blurFilter, 'same');

	% Make it binary
    BW = imbinarize(theBlurredRFmap/max(theBlurredRFmap(:)));
    BW = imclearborder(BW);
    BW = bwareafilt(BW,1);

    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});

    spatialSupportRange = spatialSupportDegs(end)-spatialSupportDegs(1);
    theCentroid(1) = spatialSupportDegs(1) + s.Centroid(1)/numel(spatialSupportDegs)*spatialSupportRange;
    theCentroid(2) = spatialSupportDegs(1) + s.Centroid(2)/numel(spatialSupportDegs)*spatialSupportRange;
end

