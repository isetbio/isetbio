function absorptionsDensity = compute(obj, oi, varargin)
% Compute the cone absorptions, possibly for multiple instances
%
% Syntax:
%   [absorptions, current, interpFilters, meanCur] = compute(obj, oi);
%   [absorptions, current, interpFilters, meanCur] = ...
%       compute(obj, oiSequence);
%
% Description:
%    Compute the response of a @cMosaic object to an optical image or an
%    optical image sequence
%
% Inputs:
%    obj                 - A @cMosaic object
%    oi                  - An optical image, or oiSequence. 
%
% Outputs:
%
% Optional key/value pairs:
%    'emPath'            - Eye movement path (Nx2 matrix). Default is [0 0].
%                          N controls the number of response instances
%                          generated.

    % Match the wavelength support of the object to that of the input optical image
    % First save the current value so we can restore it once we are done
    % with this computation
    originalValues.wave = obj.wave;
    
    % Set the new wavelength support. The cMosaic object has listeners on wave
    % that set the wavelength support of the attached pigment and macular
    % properties
    obj.wave = oiGet(oi, 'wave');
    
    % Adjust retinal irradiance to account for ecc-based variations in
    % macular pigment density
    [photons, r, c] = photonsAdjustedForEccVariationInMacularPigment(obj,oi);
    
    
    % Compute wavelength-spacing scaled quantal efficiencies
    scaledQE = obj.qe * oiGet(oi, 'bin width');
    
    % Compute density of cone absosprions, by integrating photons over
    % wavelength. This is a rectangular contiguous representation of
    % aborptions with spatial sampling equal to that of the optical image
    % The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
    absorptionsDensity = XW2RGBFormat(photons *  scaledQE, r, c);
    
    % Cones integrate photons over their aperture. This is equivalent to
    % blurring. Physically, this operation occurs as the light is gathered by the cone, 
    % before integrating over wavelength. We treat the cone aperature as
    % independent of wavelength, however, so the convolution commutes with the
    % summation over wavelength used above to get the isomerization density for
    % each class of cone. It's faster to convolve here, since there are fewer
    % bands (equal to the number of cone types) to deal with. The optical
    % image is typically sampled at 30+ wavelengths
    meanLevelForRendering = absorptionsDensity(1,1,1);
    absorptionsDensity = blurByConeAperture(obj, oi, absorptionsDensity);
    

    % All done. Restore original values
    obj.wave = originalValues.wave;

end

function  displayAbsorptionsDensity(ax, oi,absorptionsDensity, meanLevel)
    
    oiResMicrons = oiGet(oi, 'height spatial resolution')*1e6;
    oiSize = oiGet(oi, 'size');
    % Remove 10% (5% on each size to avoid edge artifacts)
    removedMargin = round(0.5*oiSize(1)*0.1);
    absorptionsDensity= absorptionsDensity(removedMargin:end-removedMargin-1, removedMargin:end-removedMargin-1,:);
    oiSize = oiSize - 2*removedMargin;
    
    oiSpatialSupportX = (1:oiSize(2));
    oiSpatialSupportX = (oiSpatialSupportX - mean(oiSpatialSupportX))*oiResMicrons;
    
    oiSpatialSupportY = (1:oiSize(1));
    oiSpatialSupportY = (oiSpatialSupportY - mean(oiSpatialSupportY))*oiResMicrons;
    
    
    contrastImage = (absorptionsDensity - meanLevel)/meanLevel;

    
    % Render
    imagesc(ax, oiSpatialSupportX, oiSpatialSupportY,contrastImage);
    set(ax, 'CLim', [-1 1], 'FontSize', 16);
    ylabel('space (microns)');
    axis(ax, 'xy');
    axis(ax, 'square');
    colormap(ax,gray(1024));
    drawnow;
end


function absorptionsDensityConvolutionImage = blurByConeAperture(obj, oi, absorptionsDensityImage)

    % Determine resolution of optical image
    oiResMicrons = oiGet(oi, 'height spatial resolution')*1e6;
    oiSize = oiGet(oi, 'size');
    
    % Compute cone apertures and their range
    coneApertureDiametersMicrons = obj.coneRFspacingsMicrons * obj.coneApertureToDiameterRatio;
    
    if (obj.eccVaryingConeBlur) && (sum(obj.eccentricityDegs) ~= 0)
        fprintf(2, 'eccVaryingConeBlur not performed for off axis cMosaics. Will blur using median cone aperture.\n');
        % Convolve using the median cone aperture across the mosaic
        absorptionsDensityConvolutionImage = ...
                convolveWithAperture(absorptionsDensityImage, median(coneApertureDiametersMicrons), oiResMicrons);
        return;
    end
    
    if (obj.eccVaryingConeBlur == false)
        % Convolve using the median cone aperture across the mosaic
        absorptionsDensityConvolutionImage = ...
                convolveWithAperture(absorptionsDensityImage, median(coneApertureDiametersMicrons), oiResMicrons);
        return;
    end
    
    tic
    % Partition optical image to subergions that correspond to different zones of cone
    % aperture size.  First define the zones by discritizing cone aperture size.
    coneApertureMicronsStepSize = 0.02;
    [coneApertureDiameterMicronsZones, ...  % the median cone aperture in this zone
        coneIndicesInZones, ...             % the IDs of cones in this zone
        nearestConeIndexTable ...           % vector storing the code ID that is nearest to each pixel of the OI
    ] = partitionOpticalImageInZones(obj.coneRFpositionsMicrons, ...
                                        coneApertureDiametersMicrons, ...
                                        coneApertureMicronsStepSize, ...
                                        oiResMicrons, oiSize);
    fprintf('Paritioning took %2.2f minutes\n', toc/60);
    
    % Allocate memory for the merged absorptions density image 
    nRows = size(absorptionsDensityImage,1);
    mCols = size(absorptionsDensityImage,2);
    absorptionsDensityMergedConvolutionsImage = bsxfun(@plus, ...
                zeros(nRows*mCols, size(absorptionsDensityImage,3)), 0*absorptionsDensityImage(1,1,:));
    

    % Perform separate convolutions with kernels based on each of
    % the discritized cone aperture zones
    for zoneIndex = 1:numel(coneApertureDiameterMicronsZones)
        fprintf('Computing blur in zone %d of %d\n', zoneIndex, numel(coneApertureDiameterMicronsZones));
        % Retrieve IDs of cones in this aperture zone
        coneIDsInZone = coneIndicesInZones{zoneIndex};
        % Retrieve median aperture of cones in this partition
        coneApertureMicronsInZone = coneApertureDiameterMicronsZones(zoneIndex);
        
        tic
        % See if we need to perform a convolution. 
        if (coneApertureDiameterMicronsZones(zoneIndex) > oiResMicrons)
            % Blur with cone aperture
            absorptionsDensityConvolutionImage = ...
                convolveWithAperture(absorptionsDensityImage, coneApertureMicronsInZone, oiResMicrons);
        else
            % Since oiResMicrons is > cone aperture, skip blurring with cone aperture
            fprintf(2,'Cone aperture is smaller that oiRes. Will not convolve with cone aperture.\n');
            absorptionsDensityConvolutionImage = absorptionsDensityImage;
        end
        fprintf('Colvolution took %2.3f seconds', toc)
        
        tic
        % Determine which pixels in the OI should receive this blur. These
        % are pixels that lie the closest than all other pixels to the cones in  this zone
        affectedOIpixels = [];
        for iCone = 1:numel(coneIDsInZone)
            affectedOIpixels = cat(2, affectedOIpixels, find(nearestConeIndexTable == coneIDsInZone(iCone)));
        end

        meanLevel = absorptionsDensityImage(1,1,1);
        % Assign these pixels from this blurred image to the final oiImage
        for coneType = 1:size(absorptionsDensityImage,3)
            channelImage = squeeze(absorptionsDensityConvolutionImage(:,:,coneType));
            absorptionsDensityMergedConvolutionsImage(affectedOIpixels, coneType) = channelImage(affectedOIpixels);
        end  
        
        fprintf('Pixel assignment took %2.2f seconds\n', toc);
        
        if (1==2)
            figure(100+zoneIndex);
            tmp = squeeze(absorptionsDensityMergedConvolutionsImage(:, 1));
            ax = subplot('Position', [0.04 0.04 0.95 0.95]);
            displayAbsorptionsDensity(ax,oi,reshape(tmp, [nRows mCols]), meanLevel);
        end
        
    end % zoneIndex
    
    % Reshape the merged convolutions image
    absorptionsDensityConvolutionImage = 0*absorptionsDensityImage;
    for coneType = 1:size(absorptionsDensityImage,3)
       absorptionsDensityVector = squeeze(absorptionsDensityMergedConvolutionsImage(:, coneType));
       absorptionsDensityConvolutionImage(:,:,coneType) = reshape(absorptionsDensityVector, [nRows mCols 1]);
    end
end

function  [coneApertureDiameterMicronsZones, coneIndicesInZones, nearestConeIndexTable] = ...
    partitionOpticalImageInZones(coneRFpositionsMicrons, coneApertureDiametersMicrons, coneApertureMicronsStepSize, oiResMicrons, oiSize)
        
    % Since optical images are always zero centered, compute center the cone
    % positions by subtracting the center of the cone mosaic
    conePosMicronsCenteredOnOpticalImage = ...
        bsxfun(@minus,coneRFpositionsMicrons, mean(coneRFpositionsMicrons,1));
    
    % Compute nearestOIpixel for each cone and 
    % nearestConeIndex for each optical image pixel
    % This operation is computationally expensive. See if we can cache the
    % computed table and reuse when the 
    %   - conePosMicronsCenteredOnOpticalImage, - oiResMicrons, - oiSize 
    % are all unchanged
    nearestConeIndexTable = ...
        opticalImagePixelsInProximityToConePositions(conePosMicronsCenteredOnOpticalImage, oiResMicrons, oiSize);
    
    % Discritize cone apertures range in zones with minimum
    % separation equal to coneApertureMicronsStepSize 
    prctileRange = prctile(coneApertureDiametersMicrons, [1 99]);
    nStepsMax = round((prctileRange(2)-prctileRange(1))/coneApertureMicronsStepSize);
    for nStepsTested = 2:nStepsMax
        coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nStepsTested);
        firstStep = coneApertureDiscritization(2)-coneApertureDiscritization(1);
        d(nStepsTested-1) = abs(firstStep-coneApertureMicronsStepSize);
    end
    [~, idx] = min(d);
    nSteps = idx+1;
   
    coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nSteps);
    fprintf('To achieve min cone aperture step of %2.3f (%d), we will discretized with %d steps\n', ...
        coneApertureMicronsStepSize, coneApertureDiscritization(2)-coneApertureDiscritization(1), nSteps);
    
    coneApertureDiscritization = cat(2, min(coneApertureDiametersMicrons), coneApertureDiscritization);
    coneApertureDiscritization = cat(2, coneApertureDiscritization, 1.01*max(coneApertureDiametersMicrons));
    
    % Determine indices of cones for each aperture zone
    coneApertureDiameterMicronsZones = zeros(1,numel(coneApertureDiscritization)-1);
    coneIndicesInZones = cell(1, numel(coneApertureDiscritization)-1);
    for i = 1:numel(coneApertureDiameterMicronsZones)
        coneApertureDiameterMicronsZones(i) = 0.5*(coneApertureDiscritization(i)+coneApertureDiscritization(i+1));
        coneIndicesInZones{i} = find(...
            (coneApertureDiametersMicrons >= coneApertureDiscritization(i)) & ...
            (coneApertureDiametersMicrons < coneApertureDiscritization(i+1)));
        coneApertureDiameterMicronsZones(i) = median(coneApertureDiametersMicrons(coneIndicesInZones{i}));
    end
end


function nearestConeIndex = opticalImagePixelsInProximityToConePositions(conePositionsMicrons, oiResMicrons, oiSize)

    oiYPosMicrons = (1:oiSize(1))*oiResMicrons;
    oiXPosMicrons = (1:oiSize(2))*oiResMicrons;
    [oiPositionsMicronsX, oiPositionsMicronsY] = ...
        meshgrid(oiXPosMicrons - mean(oiXPosMicrons), oiYPosMicrons - mean(oiYPosMicrons));
    
    oiPositionsMicrons = [oiPositionsMicronsX(:) oiPositionsMicronsY(:)];
    
    % Compute the nearest OIpixel for each cone
    %[~,nearestOIpixelIndex] = pdist2(oiPositionsMicrons, conePositionsMicrons, 'euclidean','Smallest',1);
    %assert(numel(nearestOIpixelIndex) == size(conePositionsMicrons,1), ...
    %    sprintf('Did not compute nearest OI pixel for all cones'));
    
    % Compute the nearest cone index for each OI pixel
    [~,nearestConeIndex] = pdist2(conePositionsMicrons, oiPositionsMicrons, 'euclidean','Smallest',1);
    assert(numel(nearestConeIndex) == prod(oiSize), ...
        sprintf('Did not compute nearest cone for all OI pixels'));
end

function absorptionsDensity = convolveWithAperture(absorptionsDensity, coneApertureDiameterMicrons, oiResMicrons)
    % Compute convolution kernel size in pixels
    apertureSamples = ceil(coneApertureDiameterMicrons / oiResMicrons);
    
    % Make sure it is odd
    if (mod(apertureSamples, 2) == 0)
        apertureSamples = apertureSamples + 1;
    end
    
    % Generate aperture kernel: pillnox with unit volume
    apertureKernel = zeros(apertureSamples, apertureSamples);
    [xx,yy] = sample2space(1:apertureSamples, 1:apertureSamples, oiResMicrons, oiResMicrons);
    [apertureRows, apertureCols] = meshgrid(xx,yy);
    apertureRadii = sqrt(apertureRows .^ 2 + apertureCols .^ 2);
    apertureKernel(apertureRadii <= 0.5*coneApertureDiameterMicrons) = 1;
    apertureKernel = apertureKernel ./ sum(apertureKernel(:));
    
    % Do the convolution for each cone type.
    for iCone = 1:size(absorptionsDensity,3)
        absorptionsDensity(:, :, iCone) = conv2(absorptionsDensity(:, :, iCone), apertureKernel, 'same');
    end 

end


function [photons, r, c] = photonsAdjustedForEccVariationInMacularPigment(obj,oi)
    % Retrieve retinal irradiance in photons
    photons = oiGet(oi, 'photons');
    
    % Reshape the photons for efficient computations
    [photons, r, c] = RGB2XWFormat(photons);
    
    xDegs = (0:(c-1))/c * oiGet(oi, 'h fov');
    yDegs = (0:(r-1))/r * oiGet(oi, 'v fov');
    xDegs = xDegs-mean(xDegs); yDegs = yDegs-mean(yDegs);
    [X,Y] = meshgrid(xDegs,yDegs); X = X(:); Y = Y(:); 
    eccDegs = sqrt(X.^2+Y.^2);
    
    % Extract default MP transmittance
    defaultMacularPigmentTransmittance = obj.macular.transmittance;
    
    % Compute ecc-based MP optical densities
    eccBasedMacularPigmentDensities = Macular.eccDensity(eccDegs);
    % And corresponding transmittances
    eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * obj.macular.unitDensity');

    % Compute boost factor for optical image photons so as to counteract the
    % increased transmittance through the macular pigment at increasing eccentricities 
    % due to the reduction in the MP density with eccentricity
    opticalImageBoostFactor = bsxfun(@rdivide, eccBasedMacularPigmentTransmittances, defaultMacularPigmentTransmittance');
    
    % Boost retinal image photons 
    photons = photons .* opticalImageBoostFactor;
end

