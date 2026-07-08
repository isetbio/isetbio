function computeConeAperturesRefactored(obj, lowOpticalImageResolutionWarning)
% Compute apertures from the cone positions stored in obj.
%
% ---
% Refactored for clarity by Gemini. Claims speed.  Not tested.
% ---
%
% Original logic:
% 1. Compute spacings.
% 2. Optionally smooth.
% 3. Compute ecc-varying apertures.
% 4. A deeply nested if/else block to decide whether to *use* the
%    ecc-varying apertures, a median, a min, or an imported value.
% 5. This nested block also decided on blur zones.
% 6. A final check for 'anchorAllEccVaryingParamsToTheirFovealValues'
%    could overwrite everything.
%
% New logic:
% 1. Compute spacings.
% 2. Optionally smooth.
% 3. Compute a *candidate* set of ecc-varying apertures.
% 4. (Optional) Apply rod-intrusion to this candidate set.
% 5. Use a single, flat if/elseif/else block to determine the
%    *final* apertures and blur zones based on override properties
%    (imported, anchor, eccVarying)
% 6. Assign final values and convert to degrees.
% ---
%
% To test and decide what to do.  Will need comments.
%
%-------------------------------------
% Compute apertures from the cone positions stored in obj
%
% Syntax
%  computeConeApertures(cMosaic,lowResolutionWarning)
%  cMosaic.computeConeApertures(lowResolutionWarning)
%
% Inputs:
%   lowResolutionWarning
%
% Description
%   Hairy code that NC wrote.  Will try to comment (BW). A lot of the work
%   is done in the included method
%   coneZonesFromApertureSizeAndOIresolutionve
%
% Author: NC
%
% See also
%   cMosaic
%   Includes two methods: coneZonesFromApertureSizeAndOIresolution and
%  

nConesNum = size(obj.coneRFpositionsMicrons, 1);

%% --- 1. Compute Base Cone Spacing ---
[rawConeSpacings, nearbyConeIndices] = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);

%% --- 2. (Optional) Smooth Spacings ---
if (~obj.anchorAllEccVaryingParamsToTheirFovealValues) && ...
   (isfield(obj.coneApertureModifiers, 'smoothLocalVariations')) && ...
   (obj.coneApertureModifiers.smoothLocalVariations)

    % Smooth variations in cone apertures
    smoothedSpacings = zeros(size(rawConeSpacings));
    unsmoothedRFspacingsMicrons = rawConeSpacings;
    parfor coneIndex = 1:nConesNum
        smoothedSpacings(coneIndex) = median(unsmoothedRFspacingsMicrons(nearbyConeIndices(:, coneIndex)));
    end
    obj.coneRFspacingsMicrons = smoothedSpacings;
else
    obj.coneRFspacingsMicrons = rawConeSpacings;
end

%% --- 3. Compute Initial (Candidate) Ecc-Varying Apertures ---

% Set conversion factor based on aperture shape
if (isfield(obj.coneApertureModifiers, 'shape')) && ...
   (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')) && ...
   (isfield(obj.coneApertureModifiers, 'sigma'))
    
    obj.coneApertureToDiameterRatio = 1;
    obj.coneApertureToConeCharacteristicRadiusConversionFactor = obj.coneApertureModifiers.sigma * sqrt(2.0);
else
    obj.coneApertureToConeCharacteristicRadiusConversionFactor = 0.204 * sqrt(2.0);
end

% Calculate the initial, eccentricity-varying apertures
% We treat this as a candidate value, 'eccVaryingApertures', which may or
% may not be used depending on the object's properties.
eccVaryingApertures = obj.coneRFspacingsMicrons * ...
                      obj.coneDiameterToSpacingRatio * ...
                      obj.coneApertureToDiameterRatio;

%% --- 4. (Optional) Apply Rod Intrusion Adjustment ---
% This adjustment only applies if we are in an eccentricity-varying mode.
if (obj.eccVaryingConeAperture)
    if (islogical(obj.rodIntrusionAdjustedConeAperture)) && (obj.rodIntrusionAdjustedConeAperture)
        % Compute and apply shrinkage factors
        if (isempty(obj.coneApertureRodIntrusionInducedShrinkageFactors))
            obj.computeConeApertureRodIntrusionInducedShrinkageFactors();
        end
        eccVaryingApertures = eccVaryingApertures .* obj.coneApertureRodIntrusionInducedShrinkageFactors;
        
    elseif (isnumeric(obj.rodIntrusionAdjustedConeAperture)) && (obj.rodIntrusionAdjustedConeAperture < 1)
        % Apply a simple scalar shrinkage
        eccVaryingApertures = eccVaryingApertures * obj.rodIntrusionAdjustedConeAperture;
    end
end

%% --- 5. Determine Final Apertures and Blur Zones ---
% This block replaces the complex nested if/else structure from the
% original. We check for overrides first, then handle the branching logic
% in a clear, top-down manner.

if (~isempty(obj.importedBlurDiameterMicrons))
    % --- CASE 1: Override all with imported value ---
    % All apertures and blur zones are set to the imported scalar.
    finalApertures = ones(nConesNum, 1) * obj.importedBlurDiameterMicrons;
    blurApertureDiameterMicronsZones(1) = obj.importedBlurDiameterMicrons;
    coneIndicesInZones{1} = 1:nConesNum;

elseif (obj.anchorAllEccVaryingParamsToTheirFovealValues)
    % --- CASE 2: Anchor all to foveal (min) value ---
    % All apertures and blur zones are set to the minimum (foveal) aperture.
    fovealAperture = min(eccVaryingApertures);
    finalApertures = ones(nConesNum, 1) * fovealAperture;
    blurApertureDiameterMicronsZones(1) = fovealAperture;
    coneIndicesInZones{1} = 1:nConesNum;
    fprintf('cone aperture diameter set to %f microns for all cones\n', fovealAperture);
    fprintf('cone aperture blur diameter set to %f microns for all cones\n', fovealAperture);

elseif (~obj.eccVaryingConeAperture)
    % --- CASE 3: Uniform (median) apertures ---
    % Ecc-varying is OFF. All apertures and blur zones set to the median.
    medianAperture = median(eccVaryingApertures);
    finalApertures = ones(nConesNum, 1) * medianAperture;
    blurApertureDiameterMicronsZones(1) = medianAperture;
    coneIndicesInZones{1} = 1:nConesNum;

else
    % --- CASE 4: Eccentricity-Varying Apertures ---
    % This is the default, most complex case.
    finalApertures = eccVaryingApertures;
    
    if (obj.eccVaryingConeBlur)
        % Ecc-Varying Apertures AND Ecc-Varying Blur
        % We must partition the cones into blur zones.
        if (isempty(obj.oiResMicronsForZoning))
            obj.oiResMicronsForZoning = 0.1;
        end
        
        [blurApertureDiameterMicronsZones, coneIndicesInZones] = ...
            obj.coneZonesFromApertureSizeAndOIresolution(finalApertures, lowOpticalImageResolutionWarning);
        
        if (numel(blurApertureDiameterMicronsZones) > 1)
            fprintf('Using %d blur zones\n', numel(blurApertureDiameterMicronsZones));
        end
        
    else
        % Ecc-Varying Apertures, but UNIFORM Blur
        % Apertures vary, but all are blurred with a single (median) kernel.
        blurApertureDiameterMicronsZones(1) = median(finalApertures);
        coneIndicesInZones{1} = 1:nConesNum;
    end
end

%% --- 6. Assign Final Values to Object ---
obj.coneApertureDiametersMicrons = finalApertures;
obj.blurApertureDiameterMicronsZones = blurApertureDiameterMicronsZones;
obj.coneIndicesInZones = coneIndicesInZones;

%% --- 7. Convert to Degrees ---
% Compute mean eccentricity for each zone
eccZones = zeros(1, numel(coneIndicesInZones));
for zoneIndex = 1:numel(coneIndicesInZones)
    coneIndices = coneIndicesInZones{zoneIndex};
    if isempty(coneIndices)
        eccZones(zoneIndex) = 0;
    else
        meanPosition = mean(obj.coneRFpositionsMicrons(coneIndices,:), 1);
        eccZones(zoneIndex) = sqrt(sum(meanPosition.^2));
    end
end

% Convert all micron-based values to degrees
eccRadiiMicrons = (sqrt(sum(obj.coneRFpositionsMicrons.^2, 2)))';
obj.coneRFspacingsDegs = obj.sizeMicronsToSizeDegreesForCmosaic(obj.coneRFspacingsMicrons, eccRadiiMicrons);
obj.coneApertureDiametersDegs = obj.sizeMicronsToSizeDegreesForCmosaic(obj.coneApertureDiametersMicrons, eccRadiiMicrons);
obj.blurApertureDiameterDegsZones = obj.sizeMicronsToSizeDegreesForCmosaic(obj.blurApertureDiameterMicronsZones, eccZones);

end


%%
function [coneApertureDiameterMicronsZoneBands, ...
          coneIndicesInZoneBands ...
         ] = coneZonesFromApertureSizeAndOIresolution(obj, coneApertureDiametersMicrons, lowOpticalImageResolutionWarning)

    oiResMicrons = obj.oiResMicronsForZoning;

    % Check whether there exists a cachedConePartition that is not stale
    if (obj.cachedConePartitionIsValid()) % Use refactored cache check
        fprintf('Using cached cone partition with oiRes = %f.\n', oiResMicrons);
        % Retrieve data from cache
        coneApertureDiameterMicronsZoneBands = obj.cachedConePartition.coneApertureDiameterMicronsZoneBands;
        coneIndicesInZoneBands = obj.cachedConePartition.coneIndicesInZoneBands;
        return;
    end

    fprintf('Cache is not valid. Recomputing cone aperture zones.\n');
    
    % --- 1. Discretize cone apertures into fine-grained zones ---
    coneApertureMicronsStepSize = oiResMicrons;
    prctileRange = prctile(coneApertureDiametersMicrons, [1 99]);

    if (oiResMicrons > 2 * prctileRange(1))
        fprintf(2, 'ecc-dependent blur requested but the optical image resolution is too low for this.\n');
        fprintf(2, 'Min cone aperture: %2.2f um, max: %2.2f um, oiRes: %2.2f  microns\n', ...
            prctileRange(1), prctileRange(2), oiResMicrons);
        fprintf(2, 'Suggested action: increase the number of pixels in the scene.\n');
    end

    nStepsMax = round((prctileRange(2) - prctileRange(1)) / coneApertureMicronsStepSize);
    if (nStepsMax < 2)
        nSteps = 2;
    else
        % Find nSteps for logspace such that the first step is
        % closest to coneApertureMicronsStepSize
        d = zeros(1, nStepsMax-1);
        for nStepsTested = 2:nStepsMax
            coneApertureDiscritization = logspace(log10(prctileRange(1)), log10(prctileRange(2)), nStepsTested);
            firstStep = coneApertureDiscritization(2) - coneApertureDiscritization(1);
            d(nStepsTested - 1) = abs(firstStep - coneApertureMicronsStepSize);
        end
        [~, idx] = min(d);
        nSteps = idx + 1;
    end

    coneApertureDiscritization = logspace(log10(prctileRange(1)), log10(prctileRange(2)), nSteps);
    coneApertureDiscritization = cat(2, min(coneApertureDiametersMicrons), coneApertureDiscritization);
    coneApertureDiscritization = cat(2, coneApertureDiscritization, 1.01 * max(coneApertureDiametersMicrons));

    % Determine indices of cones for each aperture zone
    zonesNum = numel(coneApertureDiscritization) - 1;
    coneApertureDiameterMicronsZones = zeros(1, zonesNum);
    coneIndicesInZones = cell(1, zonesNum);
    
    for i = 1:zonesNum
        coneIndicesInZones{i} = find(...
            (coneApertureDiametersMicrons >= coneApertureDiscritization(i)) & ...
            (coneApertureDiametersMicrons < coneApertureDiscritization(i + 1)));
        
        if isempty(coneIndicesInZones{i})
            coneApertureDiameterMicronsZones(i) = 0.5 * (coneApertureDiscritization(i) + coneApertureDiscritization(i + 1));
        else
            coneApertureDiameterMicronsZones(i) = median(coneApertureDiametersMicrons(coneIndicesInZones{i}));
        end
    end

    % --- 2. Generate aperture kernels for ALL zones ---
    % This is the expensive step.
    apertureKernelsForAllZones = cell(1, zonesNum);
    for zoneIndex = 1:zonesNum
        if isempty(coneIndicesInZones{zoneIndex})
            apertureKernelsForAllZones{zoneIndex} = [];
        else
            theKernel = obj.generateApertureKernel(coneApertureDiameterMicronsZones(zoneIndex), ...
                oiResMicrons, lowOpticalImageResolutionWarning);
            apertureKernelsForAllZones{zoneIndex} = theKernel;
        end
    end

    % --- 3. Merge zones with identical kernels into "bands" ---
    % This logic is N^2 but is much clearer than the original.
    coneIndicesInZoneBands = {};
    coneApertureDiameterMicronsZoneBands = [];
    visited = false(1, zonesNum);
    
    for zoneIndex = 1:zonesNum
        if visited(zoneIndex) || isempty(coneIndicesInZones{zoneIndex})
            continue;
        end
        
        % This is the start of a new band
        currentBandIndices = coneIndicesInZones{zoneIndex};
        apertureKernel1 = apertureKernelsForAllZones{zoneIndex};
        visited(zoneIndex) = true;
        
        % Find all other unvisited zones with an identical kernel
        for otherZoneIndex = (zoneIndex + 1):zonesNum
            if visited(otherZoneIndex) || isempty(coneIndicesInZones{otherZoneIndex})
                continue;
            end
            
            apertureKernel2 = apertureKernelsForAllZones{otherZoneIndex};
            
            % Check for identical kernels
            if isequal(size(apertureKernel1), size(apertureKernel2))
                if max(abs(apertureKernel1(:) - apertureKernel2(:))) == 0
                    % Kernels are identical, add to band
                    currentBandIndices = [currentBandIndices, coneIndicesInZones{otherZoneIndex}]; %#ok<AGROW>
                    visited(otherZoneIndex) = true;
                end
            end
        end
        
        % Save the completed band
        coneIndicesInZoneBands{end+1} = currentBandIndices; %#ok<AGROW>
        coneApertureDiameterMicronsZoneBands(end+1) = median(coneApertureDiametersMicrons(currentBandIndices)); %#ok<AGROW>
    end

    % --- 4. Save to cache ---
    obj.cachedConePartition.nConesNum = size(obj.coneRFpositionsMicrons, 1);
    obj.cachedConePartition.coneApertureModifiers = obj.coneApertureModifiers;
    obj.cachedConePartition.oiResMicronsForZoning = obj.oiResMicronsForZoning;
    obj.cachedConePartition.coneApertureDiameterMicronsZoneBands = coneApertureDiameterMicronsZoneBands;
    obj.cachedConePartition.coneIndicesInZoneBands = coneIndicesInZoneBands;
end



%% Check if the cached cone partition is valid
function isValid = cachedConePartitionIsValid(obj)

    if (isempty(obj.cachedConePartition)) || ...
       (~isfield(obj.cachedConePartition, 'oiResMicronsForZoning')) || ...
       (~isfield(obj.cachedConePartition, 'coneApertureModifiers')) || ...
       (~isfield(obj.cachedConePartition, 'nConesNum'))
        isValid = false;
        return;
    end

    % --- Perform faster, explicit checks ---
    
    % Check 1: Number of cones
    if (size(obj.coneRFpositionsMicrons, 1) ~= obj.cachedConePartition.nConesNum)
        isValid = false;
        return;
    end
    
    % Check 2: OI resolution used for zoning
    if (obj.oiResMicronsForZoning ~= obj.cachedConePartition.oiResMicronsForZoning)
        isValid = false;
        return;
    end

    % Check 3: Aperture modifiers (the slow part in the original)
    % We replace the recursive struct comparison with explicit checks
    % on the fields we *know* matter (from the main function).
    s1 = obj.coneApertureModifiers;
    s2 = obj.cachedConePartition.coneApertureModifiers;

    % Check 'shape' field
    s1_hasShape = isfield(s1, 'shape');
    s2_hasShape = isfield(s2, 'shape');
    if (s1_hasShape ~= s2_hasShape)
        isValid = false;
        return;
    elseif (s1_hasShape && ~strcmp(s1.shape, s2.shape))
        isValid = false;
        return;
    end

    % Check 'sigma' field
    s1_hasSigma = isfield(s1, 'sigma');
    s2_hasSigma = isfield(s2, 'sigma');
    if (s1_hasSigma ~= s2_hasSigma)
        isValid = false;
        return;
    elseif (s1_hasSigma && (s1.sigma ~= s2.sigma))
        isValid = false;
        return;
    end
    
    % Check 'smoothLocalVariations' field
    s1_hasSmooth = isfield(s1, 'smoothLocalVariations');
    s2_hasSmooth = isfield(s2, 'smoothLocalVariations');
    if (s1_hasSmooth ~= s2_hasSmooth)
        isValid = false;
        return;
    elseif (s1_hasSmooth && (s1.smoothLocalVariations ~= s2.smoothLocalVariations))
         isValid = false;
        return;
    end

    % If all checks passed
    isValid = true;
end