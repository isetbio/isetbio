function computeConeApertures(obj, lowOpticalImageResolutionWarning)

    % Compute unsmoothed spacings from positions
    [obj.coneRFspacingsMicrons, nearbyConeIndices] = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);
    
    % Retrieve number of cones
    nConesNum = size(obj.coneRFpositionsMicrons,1);
        

    if (~obj.anchorAllEccVaryingParamsToTheirFovealValues)
        if (isfield(obj.coneApertureModifiers, 'smoothLocalVariations')) && (obj.coneApertureModifiers.smoothLocalVariations)
            % Smooth variations in cone apertures
            coneRFspacingsMicronsTmp = 0*obj.coneRFspacingsMicrons;
            unsmoothedRFspacingsMicrons = obj.coneRFspacingsMicrons;
            parfor coneIndex = 1:nConesNum
                coneRFspacingsMicronsTmp(coneIndex) = median(unsmoothedRFspacingsMicrons(nearbyConeIndices(:,coneIndex)));
            end
            obj.coneRFspacingsMicrons = coneRFspacingsMicronsTmp;
            clear 'coneRFspacingsMicronsTmp'
        else
            %fprintf(2, 'Will NOT smooth local variations in cone aperture/spacing\n');
        end
    end

    
    if (...
            (isfield(obj.coneApertureModifiers, 'shape')) && ...
            (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')) && ...
            (isfield(obj.coneApertureModifiers, 'sigma')) ...
        )
        obj.coneApertureToDiameterRatio = 1;
        obj.coneApertureToConeCharacteristicRadiusConversionFactor = obj.coneApertureModifiers.sigma * sqrt(2.0);
    else
        obj.coneApertureToConeCharacteristicRadiusConversionFactor = 0.204*sqrt(2.0);
    end
    
    % Compute cone aperture diameters based on their local spacing and the coneDiameterToSpacingRatio
    coneAperturesMicrons = obj.coneRFspacingsMicrons * ...
                           obj.coneDiameterToSpacingRatio * ... 
                           obj.coneApertureToDiameterRatio;

    reportZonesUsed = false;
  
    if (obj.eccVaryingConeAperture)

        % Rod intrusion adjustment
        if (islogical(obj.rodIntrusionAdjustedConeAperture))
            if (obj.rodIntrusionAdjustedConeAperture)
                if (isempty(obj.coneApertureRodIntrusionInducedShrinkageFactors))
                    % Compute coneApertureRodIntrusionInducedShrinkageFactors to
                    % account for progressive rod intrusion with eccentricity which
                    % shrinks the apertures of cones (relative to cone spacing)
                    obj.computeConeApertureRodIntrusionInducedShrinkageFactors();
                end
                coneAperturesMicrons = coneAperturesMicrons .* obj.coneApertureRodIntrusionInducedShrinkageFactors;
            end
        else
            if (obj.rodIntrusionAdjustedConeAperture < 1)
                coneAperturesMicrons = coneAperturesMicrons * obj.rodIntrusionAdjustedConeAperture;
            end
        end

        if (obj.eccVaryingConeBlur)
            % oiResMicrons: The lower this is the more finely
            % cone apertures will be zoned, but compute time will increase
            if (isempty(obj.oiResMicronsForZoning))
                obj.oiResMicronsForZoning = 0.1;
            else
                reportZonesUsed = true;
            end
            %tStart = clock;
            % Kernel aperture blur - variable with eccentricity
            % Partition cones into zones based on their aperture size.
            [blurApertureDiameterMicronsZones, ...  % the median cone aperture in this zone band
            coneIndicesInZones  ...                 % the IDs of cones in this zone band
            ] = coneZonesFromApertureSizeAndOIresolution(obj, coneAperturesMicrons, lowOpticalImageResolutionWarning);
            %fprintf('Cone zoning based on aperture size took %f seconds.\n', etime(clock, tStart));
            for zoneIndex = 1:numel(coneIndicesInZones)
                eccZones(zoneIndex) = sqrt(sum((mean(obj.coneRFpositionsMicrons(coneIndicesInZones{zoneIndex},:), 1)).^2));
            end
        else
            % Kernel aperture blur - single aperture for all eccentricities
            if (isempty(obj.importedBlurDiameterMicrons))
                blurApertureDiameterMicronsZones(1) = median(coneAperturesMicrons);
            else
               % Override both blur diameter & summation aperture with imported value
               blurApertureDiameterMicronsZones(1) = obj.importedBlurDiameterMicrons;
               coneAperturesMicrons = coneAperturesMicrons * 0 + blurApertureDiameterMicronsZones(1);
            end
            coneIndicesInZones{1} = 1:nConesNum; 
            eccZones(1) = sqrt(sum((mean(obj.coneRFpositionsMicrons(coneIndicesInZones{1},:), 1)).^2));
        end
    else
        % Kernel aperture blur - single aperture for all eccentricities
        if (isempty(obj.importedBlurDiameterMicrons))
            if (obj.anchorAllEccVaryingParamsToTheirFovealValues)
                blurApertureDiameterMicronsZones(1) = min(coneAperturesMicrons);
            else
                blurApertureDiameterMicronsZones(1) = median(coneAperturesMicrons);
            end
        else
            % Override both blur diameter & summation aperture with imported value
            blurApertureDiameterMicronsZones(1) = obj.importedBlurDiameterMicrons;
        end
        coneAperturesMicrons = coneAperturesMicrons*0 + blurApertureDiameterMicronsZones(1);
        coneIndicesInZones{1} = 1:nConesNum;
        eccZones(1) = sqrt(sum((mean(obj.coneRFpositionsMicrons(coneIndicesInZones{1},:), 1)).^2));
    end
    
    if (reportZonesUsed)
        fprintf('Using %d blur zones\n', numel(blurApertureDiameterMicronsZones));
    end
    
    % Attach to cMosaic object
    if (obj.anchorAllEccVaryingParamsToTheirFovealValues)
        obj.coneApertureDiametersMicrons = coneAperturesMicrons*0 + min(coneAperturesMicrons);
        obj.blurApertureDiameterMicronsZones = blurApertureDiameterMicronsZones*0 + min(coneAperturesMicrons);
        fprintf('cone aperture diameter set to %f microns for all cones\n', min(coneAperturesMicrons));
        fprintf('cone aperture blur diameter set to %f microns for all cones\n', min(blurApertureDiameterMicronsZones));
        obj.coneIndicesInZones{1} = 1:nConesNum;
    else
        obj.coneApertureDiametersMicrons = coneAperturesMicrons;
        obj.blurApertureDiameterMicronsZones = blurApertureDiameterMicronsZones;
        obj.coneIndicesInZones = coneIndicesInZones;
    end

    
    % Convert to degs
    eccRadiiMicrons = (sqrt(sum(obj.coneRFpositionsMicrons.^2,2)))';
    obj.coneRFspacingsDegs = obj.sizeMicronsToSizeDegreesForCmosaic(obj.coneRFspacingsMicrons, eccRadiiMicrons);
    obj.coneApertureDiametersDegs = obj.sizeMicronsToSizeDegreesForCmosaic(obj.coneApertureDiametersMicrons, eccRadiiMicrons);
    obj.blurApertureDiameterDegsZones = obj.sizeMicronsToSizeDegreesForCmosaic(obj.blurApertureDiameterMicronsZones, eccZones);
end


% Method to partition cones into zones based on cone
% aperture size and current optical image resolution
function [coneApertureDiameterMicronsZoneBands, ...     % the median cone aperture in this zone band
          coneIndicesInZoneBands  ...                   % the IDs of cones in this zone band
         ] = coneZonesFromApertureSizeAndOIresolution(obj, coneApertureDiametersMicrons, lowOpticalImageResolutionWarning)
        

    oiResMicrons = obj.oiResMicronsForZoning;
    
    % Check whether there exists a cachedConePartition that is not stale
    if (cachedConePartitionIsValid(obj))
        fprintf('Using cached cone partition with oiRes = %f.\n', oiResMicrons);
        % Retrieve data from cache
        coneApertureDiameterMicronsZoneBands = obj.cachedConePartition.coneApertureDiameterMicronsZoneBands;
        coneIndicesInZoneBands = obj.cachedConePartition.coneIndicesInZoneBands;
        return;
    end
    
    fprintf('Cache is not valid. Recomputing cone aperture zones.\n');
    
    
    % Discritize cone apertures range in zones with minimum
    % separation equal to coneApertureMicronsStepSize 
    coneApertureMicronsStepSize = oiResMicrons;
    
    prctileRange = prctile(coneApertureDiametersMicrons, [1 99]);
    
    if (oiResMicrons > 2*prctileRange(1))
        fprintf(2,'ecc-dependent blur requested but the optical image resolution is too low for this.');
        fprintf(2,'Min cone aperture: %2.2f um, max: %2.2f um, oiRes: %2.2f  microns\n', ...
            prctileRange(1), prctileRange(2), oiResMicrons);
        fprintf(2,'Suggested action: increase the number of pixels in the scene.\n');
    end

    
    nStepsMax = round((prctileRange(2)-prctileRange(1))/coneApertureMicronsStepSize);
    if (nStepsMax < 2)
        nSteps = 2;
    else
        for nStepsTested = 2:nStepsMax
            coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nStepsTested);
            firstStep = coneApertureDiscritization(2)-coneApertureDiscritization(1);
            d(nStepsTested-1) = abs(firstStep-coneApertureMicronsStepSize);
        end
        
        [~, idx] = min(d);
        nSteps = idx+1;
    end

    coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nSteps);
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
    
    % Compute aperture kernel for all zones
    zonesNum = numel(coneApertureDiameterMicronsZones);
    apertureKernelsForAllZones = cell(1, zonesNum);
    for zoneIndex = 1:zonesNum
        theKernel = obj.generateApertureKernel(coneApertureDiameterMicronsZones(zoneIndex), ...
            oiResMicrons, lowOpticalImageResolutionWarning);
        apertureKernelsForAllZones{zoneIndex} = theKernel;
    end

    % See if some aperture kernels are identical in neighboring zones and
    % if they are, merge the zones into zone bands
    mergedZones = cell(1,zonesNum);
    lastMergedZone = 0;
    for zoneIndex1 = 1:zonesNum
        apertureKernel1 = apertureKernelsForAllZones{zoneIndex1};
        startingZoneIndex2 = max([zoneIndex1, lastMergedZone]) + 1;
        for zoneIndex2 = startingZoneIndex2:zonesNum
            apertureKernel2 = apertureKernelsForAllZones{zoneIndex2};
            if (size(apertureKernel1,1) == size(apertureKernel2,1))
                kernelDiff = (apertureKernel1-apertureKernel2)/max(apertureKernel1(:))*100;
                maxDiff = max(abs(kernelDiff(:)));
                if (maxDiff == 0)
                    if (isempty(min(mergedZones{zoneIndex1})))
                        mergedZones{zoneIndex1} = [zoneIndex1 zoneIndex2];
                    else
                        mergedZones{min(mergedZones{zoneIndex1})} = cat(2, mergedZones{min(mergedZones{zoneIndex1})}, zoneIndex2);
                    end
                    lastMergedZone = zoneIndex2;
                end
            end
        end
    end

    zoneBands = {};
    lastZoneMerged = 0;
    for i = 1:zonesNum
        if (~isempty(mergedZones{i}))
            zoneBands{numel(zoneBands)+1} = mergedZones{i};
            lastZoneMerged = max(mergedZones{i});
        else
            if (i > lastZoneMerged)
                zoneBands{numel(zoneBands)+1} = i;
                lastZoneMerged = i;
            end
        end
    end
    
    % Assemble zone bands
    coneIndicesInZoneBands = cell(1, numel(zoneBands));
    coneApertureDiameterMicronsZoneBands = zeros(1, numel(zoneBands));

    zonedConeIndices = [];
    for zoneBandIndex = 1:numel(zoneBands)
        mergedZones = zoneBands{zoneBandIndex};
        coneIndices = [];
        for k = 1:numel(mergedZones)
            zoneIndex = mergedZones(k);
            coneIndices = cat(2, coneIndices, coneIndicesInZones{zoneIndex});
            zonedConeIndices = cat(2, zonedConeIndices, coneIndicesInZones{zoneIndex});
        end
        coneIndicesInZoneBands{zoneBandIndex} = coneIndices;
        coneApertureDiameterMicronsZoneBands(zoneBandIndex) = median(coneApertureDiametersMicrons(coneIndices));
    end

    % Save to cache: params that, if different, will trigger are-partitioning
    obj.cachedConePartition.nConesNum = size(obj.coneRFpositionsMicrons,1);
    obj.cachedConePartition.coneApertureModifiers = obj.coneApertureModifiers;
    obj.cachedConePartition.oiResMicronsForZoning = obj.oiResMicronsForZoning;
    
    % The partitioning data
    obj.cachedConePartition.coneApertureDiameterMicronsZoneBands = coneApertureDiameterMicronsZoneBands;
    obj.cachedConePartition.coneIndicesInZoneBands = coneIndicesInZoneBands;
    
end


function isValid = cachedConePartitionIsValid(obj)
    
    
    if (isempty(obj.cachedConePartition)) || ...
       (~isfield(obj.cachedConePartition, 'oiResMicronsForZoning')) || ...
       (~isfield(obj.cachedConePartition, 'coneApertureModifiers')) || ...
       (~isfield(obj.cachedConePartition, 'nConesNum')) 
        isValid = false;
    else
        % Check for properties that may affect the partitioning of cone
        % apertures into zones
        isValid = (...
            (size(obj.coneRFpositionsMicrons,1) == obj.cachedConePartition.nConesNum) && ...
            (isempty(RecursivelyCompareStructs('s1', obj.coneApertureModifiers, 's2', obj.cachedConePartition.coneApertureModifiers, ...
                  'failOnMissingFields', true, 'graphMismatchedData', false, 'verbosityLevel', 0))) && ...
            (obj.oiResMicronsForZoning == obj.cachedConePartition.oiResMicronsForZoning) ...
            );
    end
    
end

