% Method to partition cones into zones based on cone
% aperture size and current optical image resolution
function [coneApertureDiameterMicronsZoneBands, ...     % the median cone aperture in this zone band
          coneIndicesInZoneBands  ...                   % the IDs of cones in this zone band
            ] = coneZonesFromApertureSizeAndOIresolution(obj, coneApertureDiametersMicrons, oiResMicrons)
        
        
    % Check whether there exists a cachedConePartition that is not stale
    if (cachedConePartitionIsValid(obj, oiResMicrons))
        fprintf('No change in oiResMicrons, so will use cached cone partition.\n');
        % Retrieve data from cache
        coneApertureDiameterMicronsZoneBands = obj.cachedConePartition.coneApertureDiameterMicronsZoneBands;
        coneIndicesInZoneBands = obj.cachedConePartition.coneIndicesInZoneBands;
        return;
    end
    
    
    % Discritize cone apertures range in zones with minimum
    % separation equal to coneApertureMicronsStepSize 
    prctileRange = prctile(coneApertureDiametersMicrons, [1 99]);
    
    if (oiResMicrons > 2*prctileRange(1))
        fprintf(2,'ecc-dependent blur requested but the optical image resolution is too low for this.');
        fprintf(2,'Min cone aperture: %2.2f um, max: %2.2f um, oiRes: %2.2f  microns\n', ...
            prctileRange(1), prctileRange(2), oiResMicrons);
        fprintf(2,'Suggested action: increase the number of pixels in the scene.\n');
    end

    coneApertureMicronsStepSize = 0.02;
    nStepsMax = round((prctileRange(2)-prctileRange(1))/coneApertureMicronsStepSize);
    for nStepsTested = 2:nStepsMax
        coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nStepsTested);
        firstStep = coneApertureDiscritization(2)-coneApertureDiscritization(1);
        d(nStepsTested-1) = abs(firstStep-coneApertureMicronsStepSize);
    end
    [~, idx] = min(d);
    nSteps = idx+1;
    

    coneApertureDiscritization = logspace(log10(prctileRange(1)),log10(prctileRange(2)), nSteps);
    %fprintf('To achieve min cone aperture step of %2.3f (%d), we will discretized with %d steps\n', ...
    %    coneApertureMicronsStepSize, coneApertureDiscritization(2)-coneApertureDiscritization(1), nSteps);
    
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
        theKernel = cMosaic.generateApertureKernel(coneApertureDiameterMicronsZones(zoneIndex), oiResMicrons);
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

    % Save to cache
    obj.cachedConePartition.oiResMicrons = oiResMicrons;
    obj.cachedConePartition.coneApertureDiameterMicronsZoneBands = coneApertureDiameterMicronsZoneBands;
    obj.cachedConePartition.coneIndicesInZoneBands = coneIndicesInZoneBands;
    
end


function isValid = cachedConePartitionIsValid(obj, oiResMicrons)
    if (isempty(obj.cachedConePartition))
        isValid = false;
    else
        isValid = ...
            (oiResMicrons == obj.cachedConePartition.oiResMicrons);
            
    end
    
end

