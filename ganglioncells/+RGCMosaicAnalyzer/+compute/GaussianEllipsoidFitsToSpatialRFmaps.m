function theFittedGaussianEllipsoids = GaussianEllipsoidFitsToSpatialRFmaps(theRFmaps, spatialSupportDegs, ...
    GaussianEllipsoidSamplesNum, onlyFitPositiveRFmap, fixCentroid, smoothingKernelSizePixels, visualizeSmoothingAndFitting)

    fprintf('Will fit Gaussian ellispoids to the %d RF maps.\n', numel(theRFmaps));

	theFittedGaussianEllipsoids = cell(1, numel(theRFmaps));

    if (~isempty(smoothingKernelSizePixels))
        kernelSigmaPixels = smoothingKernelSizePixels/6;
        kernelSizePixels = ceil(smoothingKernelSizePixels);
        % Make it odd so the Gaussian peaks at a single pixel
        if (mod(kernelSizePixels,2) == 0)
           kernelSizePixels = kernelSizePixels + 1; 
        end 
        smoothingKernel =  fspecial('gaussian', kernelSizePixels, kernelSigmaPixels);
    else
        smoothingKernel = [];
    end


    if (visualizeSmoothingAndFitting)
        % Serial for execution since we are visualizing
        for iRGC = 1:numel(theRFmaps)
        	fprintf('Fitting Gaussian ellispoid to RF map %d of %d.\n', iRGC, numel(theRFmaps));
            theFittedGaussianEllipsoids{iRGC} = doTheFit(theRFmaps{iRGC}, spatialSupportDegs, ...
                GaussianEllipsoidSamplesNum, onlyFitPositiveRFmap, fixCentroid, ...
                smoothingKernel, visualizeSmoothingAndFitting);
        end
    else
        % Parfor execution since we are not visualizing
        parfor iRGC = 1:numel(theRFmaps)
            fprintf('Fitting Gaussian ellispoid to RF map %d of %d.\n', iRGC, numel(theRFmaps));
            theFittedGaussianEllipsoids{iRGC} = doTheFit(theRFmaps{iRGC}, spatialSupportDegs, ...
                GaussianEllipsoidSamplesNum, onlyFitPositiveRFmap, fixCentroid, ...
                smoothingKernel, visualizeSmoothingAndFitting);
        end
    end
end

function theFitStruct = doTheFit(theRFmap, spatialSupportDegs, ...
    GaussianEllipsoidSamplesNum, onlyFitPositiveRFmap, fixCentroid, ...
    smoothingKernel, visualizeSmoothingAndFitting)

    if (onlyFitPositiveRFmap)
       theRFmap(theRFmap<0) = 0;
    end

    if (visualizeSmoothingAndFitting)
            hFig = figure(33); clf;
            set(hFig, 'Position', [10 10 1000 1000]);
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);

            imagesc(ax1,spatialSupportDegs, spatialSupportDegs, theRFmap);
            axis(ax1, 'image')
            axis(ax1, 'xy');
            title(ax1, 'original RFmap (not smoothed)')
    end

    if (~isempty(smoothingKernel))
        theSmoothedRFmap = conv2(theRFmap, smoothingKernel, 'same');

        if (visualizeSmoothingAndFitting)
            imagesc(ax2,smoothingKernel);
            axis(ax2, 'image');
            axis(ax2, 'xy');
            title(ax2, 'smoothing kernel')

            imagesc(ax3, spatialSupportDegs, spatialSupportDegs, theSmoothedRFmap);
            axis(ax3,'image');
            axis(ax3, 'xy');
            title(ax3, 'fitted RFmap (smoothed)')
        end
        
        theRFmap = theSmoothedRFmap;
    end

    rfMax = max(abs(theRFmap(:)));

    assert(...
        (size(theRFmap,1) == numel(spatialSupportDegs)) && (size(theRFmap,2) == numel(spatialSupportDegs)), ...
        sprintf('size error! RFmap is %d x %d, spatialSupport is %d', size(theRFmap,1), size(theRFmap,2), numel(spatialSupportDegs)));
        

    % Estimate the RF centroid
    theCentroid = rfCentroid(theRFmap, spatialSupportDegs);

    if (fixCentroid) && (~isempty(theCentroid))
        [theFittedGaussianEllipsoidParams, ...
        theFittedGaussianEllipsoid, ...
        theFittedGaussianEllipsoidSpatialSupportXDegs, ...
        theFittedGaussianEllipsoidSpatialSupportYDegs] = RGCMosaicConstructor.helper.fit.GaussianEllipsoidTo2DRFmap(...
            spatialSupportDegs, spatialSupportDegs, theRFmap/rfMax, GaussianEllipsoidSamplesNum, ...
            'fixedCentroid', theCentroid);
    else
        [theFittedGaussianEllipsoidParams, ...
        theFittedGaussianEllipsoid, ...
        theFittedGaussianEllipsoidSpatialSupportXDegs, ...
        theFittedGaussianEllipsoidSpatialSupportYDegs] = RGCMosaicConstructor.helper.fit.GaussianEllipsoidTo2DRFmap(...
            spatialSupportDegs, spatialSupportDegs, theRFmap/rfMax, GaussianEllipsoidSamplesNum, ...
            'fixedCentroid', []);
    end

    if (visualizeSmoothingAndFitting)
        figure(hFig);

        %if (~isempty(smoothingKernel))
        
            hold(ax2, 'on')
            contour(ax2,theFittedGaussianEllipsoidSpatialSupportXDegs, theFittedGaussianEllipsoidSpatialSupportYDegs, ...
                theFittedGaussianEllipsoid*rfMax, [exp(-0.502) exp(-0.50)]*rfMax, ...
                'Color', [0 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');

 
            hold(ax3, 'on');
            contour(ax3, theFittedGaussianEllipsoidSpatialSupportXDegs, theFittedGaussianEllipsoidSpatialSupportYDegs, ...
                theFittedGaussianEllipsoid*rfMax, [exp(-0.502) exp(-0.50)]*rfMax, ...
                'Color', [0 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');
        
            ax4 = subplot(2,2,4);
            imagesc(ax4, theFittedGaussianEllipsoidSpatialSupportXDegs, theFittedGaussianEllipsoidSpatialSupportYDegs, ...
                theFittedGaussianEllipsoid);
            axis(ax4, 'image'); axis(ax4, 'xy');

            hold(ax4, 'on');
            contour(ax4, theFittedGaussianEllipsoidSpatialSupportXDegs, theFittedGaussianEllipsoidSpatialSupportYDegs, ...
                theFittedGaussianEllipsoid, [exp(-0.502) exp(-0.50)], ...
                'Color', [0 0 0], 'FaceColor', 'none', 'LineWidth', 1.5, 'LineStyle', '-');
            axis(ax4, 'image');
            drawnow;
            pause
    end

    minorSigmaDegs = min([theFittedGaussianEllipsoidParams.xSigma theFittedGaussianEllipsoidParams.ySigma]);
    majorSigmaDegs = max([theFittedGaussianEllipsoidParams.xSigma theFittedGaussianEllipsoidParams.ySigma]);

    theFitStruct = struct(...
            'rfMap', theFittedGaussianEllipsoid * rfMax, ...
            'xSupportDegs', theFittedGaussianEllipsoidSpatialSupportXDegs, ...
            'ySupportDegs', theFittedGaussianEllipsoidSpatialSupportYDegs, ...
            'minorSigmaDegs', minorSigmaDegs, ...
            'majorSigmaDegs', majorSigmaDegs, ...
            'x0', theFittedGaussianEllipsoidParams.x0, ...
            'y0', theFittedGaussianEllipsoidParams.y0, ...
            'rotationDegs', theFittedGaussianEllipsoidParams.rotationDegs ...
        );

end

function theCentroid = rfCentroid(theRFmap, spatialSupportDegs)
	blurFilter = ones(10,10); %fspecial('gaussian', 20, 4);
	theBlurredRFmap = conv2(theRFmap, blurFilter, 'same');

	% Make it binary
    BW = imbinarize(theBlurredRFmap/max(theBlurredRFmap(:)));
    %BW = imclearborder(BW);
    BW = bwareafilt(BW,5);

    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    
    if (isempty(s))
        theCentroid = [];
        return;
    end

    theCentroid = [0 0];
    for i = 1:numel(s)
        theCentroid = theCentroid + s(i).Centroid;
    end

    theCentroid = theCentroid/numel(s);
    
    spatialSupportRange = spatialSupportDegs(end)-spatialSupportDegs(1);
    theCentroid(1) = spatialSupportDegs(1) + theCentroid(1)/numel(spatialSupportDegs)*spatialSupportRange;
    theCentroid(2) = spatialSupportDegs(1) + theCentroid(2)/numel(spatialSupportDegs)*spatialSupportRange;
end

