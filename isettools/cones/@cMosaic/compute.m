function [noiseFreeAbsorptionsCount, noisyAbsorptionInstances, photoCurrents, photoCurrentInstances, responseTemporalSupport] = compute(obj, oi, varargin)
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
%    'emPaths'           - Eye movement paths. Either empty or a matrix of [nTrials x N x 2]. Default is [].
%                          

    % Parse input
    p = inputParser;
    p.addParameter('withFixationalEyeMovements', false, @islogical);
    p.addParameter('opticalImagePositionDegs', 'mosaic-centered', @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
    p.addParameter('nTimePoints', [], @isscalar);
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('seed', 1, @isnumeric);
    p.addParameter('verbosityLevel', 'none', @(x)ismember(x, {'default', 'min', 'max'}));
    p.parse(varargin{:});
    

    if (isempty(obj.minRFpositionMicrons))
        noiseFreeAbsorptionsCount = [];
        noisyAbsorptionInstances = []; 
        photoCurrents = []; 
        photoCurrentInstances = [];
        responseTemporalSupport = [];
        return;
    end
     
    % Parse optional inputs
    nTimePoints = p.Results.nTimePoints;
    nTrials = p.Results.nTrials;
    noiseSeed = p.Results.seed;
    verbosityLevel = p.Results.verbosityLevel;
    
    % Validate and decode fixational eye movement optional input
    [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
        validateAndDecodeFixationalEyeMovements(obj, p.Results.withFixationalEyeMovements, nTrials, nTimePoints);
    
    % Decode and decode opticalImagePositionDegs optional input
    opticalImagePositionMicrons = validateAndDecodeOpticalImagePosition(obj,p.Results.opticalImagePositionDegs);
    
    % Match the wavelength support of the object to that of the input optical image
    % First save the current value so we can restore it once we are done
    % with this computation
    originalValues.wave = obj.wave;
    
    % Set the new wavelength support. The cMosaic object has listeners on wave
    % that set the wavelength support of the attached pigment and macular
    % properties
    obj.wave = oiGet(oi, 'wave');
    
    % Compute wavelength-spacing scaled quantal efficiencies
    scaledQE = obj.qe * oiGet(oi, 'bin width');
    
    % Retrieve oiRes
    oiResMicrons = oiGet(oi, 'height spatial resolution')*1e6;
    
    % Retrieve wavelength support
    oiWave = oiGet(oi, 'wave');
    
    % Generate oiPositions
    spatialSupportMeters = oiGet(oi, 'spatial support');
    spatialSupportXMicrons = squeeze(spatialSupportMeters(1,1:end,1)) * 1e6;
    spatialSupportYMicrons = squeeze(spatialSupportMeters(1:end,1,2)) * 1e6;
    spatialSupportXMicrons = spatialSupportXMicrons - mean(spatialSupportXMicrons);
    spatialSupportYMicrons = spatialSupportYMicrons - mean(spatialSupportYMicrons);
    
    % Translate oi spatial support 
    spatialSupportXMicrons = spatialSupportXMicrons + opticalImagePositionMicrons(1);
    spatialSupportYMicrons = spatialSupportYMicrons + opticalImagePositionMicrons(2);
    
    % Determine if spatial support needs to be expanded so that the optical image
    % extends over the entire cone mosaic
    dx = spatialSupportXMicrons(2)-spatialSupportXMicrons(1);
    dy = spatialSupportYMicrons(2)-spatialSupportYMicrons(1);
    minEMpos = squeeze(min(emPathsMicrons,[],1));
    maxEMpos = squeeze(max(emPathsMicrons,[],1));
    
    if (obj.minRFpositionMicrons(1)+minEMpos(1) < min(spatialSupportXMicrons))
        additionalPixelsXLeft = round((min(spatialSupportXMicrons) - (obj.minRFpositionMicrons(1)+minEMpos(1)))/dx);
    else
        additionalPixelsXLeft = 0;
    end
    
    if (obj.minRFpositionMicrons(2)+minEMpos(2) < min(spatialSupportYMicrons))
        additionalPixelsYBottom = round((min(spatialSupportYMicrons) - (obj.minRFpositionMicrons(2)+minEMpos(2)))/dy);
    else
        additionalPixelsYBottom = 0;
    end
    
    if (obj.maxRFpositionMicrons(1)+maxEMpos(1) > max(spatialSupportXMicrons))
        additionalPixelsXRight = round((obj.maxRFpositionMicrons(1)+maxEMpos(1) - max(spatialSupportXMicrons))/dx);
    else
        additionalPixelsXRight = 0;
    end
    
    if (obj.maxRFpositionMicrons(2)+maxEMpos(2) > max(spatialSupportYMicrons))
        additionalPixelsYTop = round((obj.maxRFpositionMicrons(2)+maxEMpos(2) - max(spatialSupportYMicrons))/dy);
    else
        additionalPixelsYTop = 0;
    end
    
    % Extend spatial support vectors as needed
    if (additionalPixelsXLeft > 0)
        extraXPos = -fliplr((1:additionalPixelsXLeft)*dx);
        spatialSupportXMicrons = cat(2, extraXPos+spatialSupportXMicrons(1), spatialSupportXMicrons); 
    end
    if (additionalPixelsXRight > 0)
        extraXPos = (1:additionalPixelsXRight)*dx;
        spatialSupportXMicrons = cat(2, spatialSupportXMicrons, extraXPos+spatialSupportXMicrons(end)); 
    end
    
    if (additionalPixelsYBottom > 0)
        extraYPos = -fliplr((1:additionalPixelsYBottom)*dy);
        spatialSupportYMicrons = cat(1, extraYPos'+spatialSupportYMicrons(1), spatialSupportYMicrons);
    end
    if (additionalPixelsYTop > 0)
        extraYPos = (1:additionalPixelsYTop)*dy;
        spatialSupportYMicrons = cat(1,  spatialSupportYMicrons, extraYPos'+spatialSupportYMicrons(end));
    end
    
    
    if (~isempty(obj.micronsPerDegreeApproximation))
        spatialSupportXDegrees = spatialSupportXMicrons/obj.micronsPerDegreeApproximation;  
        spatialSupportYDegrees = spatialSupportYMicrons/obj.micronsPerDegreeApproximation;
    else
        spatialSupportXDegrees = RGCmodels.Watson.convert.rhoMMsToDegs(spatialSupportXMicrons*1e-3);
        spatialSupportYDegrees = RGCmodels.Watson.convert.rhoMMsToDegs(spatialSupportYMicrons*1e-3);
    end
    
    
    [oiPositionsDegsXgrid, oiPositionsDegsYgrid] = meshgrid(spatialSupportXDegrees, spatialSupportYDegrees);
    oiPositionsDegs = [oiPositionsDegsXgrid(:), oiPositionsDegsYgrid(:)];
    oiPositionsVectorsMicrons = {spatialSupportYMicrons(:), spatialSupportXMicrons(:)};
    
    % Compute cone aperture diameters based on their local spacing
    coneApertureDiametersMicrons = obj.coneRFspacingsMicrons * obj.coneApertureToDiameterRatio;
    conesNum = numel(coneApertureDiametersMicrons);
    
    if (obj.eccVaryingConeAperture)
        if (obj.eccVaryingConeBlur)
            %tStart = clock;
            % Kernel aperture blur - variable with eccentricity
            % Partition cones into zones based on their aperture size.
            [blurApertureDiameterMicronsZones, ...  % the median cone aperture in this zone band
            coneIndicesInZones  ...                 % the IDs of cones in this zone band
            ] = obj.coneZonesFromApertureSizeAndOIresolution(coneApertureDiametersMicrons, oiResMicrons);
            %fprintf('Cone zoning based on aperture size took %f seconds.\n', etime(clock, tStart));
        else
            % Kernel aperture blur - single aperture for all eccentricities
            if (isempty(obj.importedBlurDiameterMicrons))
                blurApertureDiameterMicronsZones(1) = median(coneApertureDiametersMicrons);
            else
                blurApertureDiameterMicronsZones(1) = obj.importedBlurDiameterMicrons;
            end
            coneIndicesInZones{1} = 1:conesNum; 
        end
    else
        % Kernel aperture blur - single aperture for all eccentricities
        if (isempty(obj.importedBlurDiameterMicrons))
            blurApertureDiameterMicronsZones(1) = median(coneApertureDiametersMicrons);
        else
            blurApertureDiameterMicronsZones(1) = obj.importedBlurDiameterMicrons;
        end
        coneApertureDiametersMicrons = coneApertureDiametersMicrons*0 + blurApertureDiameterMicronsZones(1);
        coneIndicesInZones{1} = 1:conesNum;
    end
    
    % Set the object's blurApertureDiameterMicronsZones property
    obj.blurApertureDiameterMicronsZones = blurApertureDiameterMicronsZones;
    
    % Retrieve retinal irradiance in photons
    photons = oiGet(oi, 'photons');
   
    % Flip the optical image upside-down because the y-coords in the
    % oi spatial support vectors increase from top -> bottom (y-coords in an image)
    % whereas cone y-positions increase from bottom -> top 
    for k = 1:size(photons,3)
        photons(:,:,k) = flipud(squeeze(photons(:,:,k)));
    end
    
    % Extend photons field as needed
    if ((additionalPixelsXLeft > 0) || (additionalPixelsXRight > 0) || (additionalPixelsYBottom > 0) || (additionalPixelsYTop > 0))
        originalPhotons = photons;
        originalXpixelsNum = size(photons,2);
        originalYpixelsNum = size(photons,1);
        newXpixelsNum = originalXpixelsNum + additionalPixelsXLeft + additionalPixelsXRight;
        newYpixelsNum = originalYpixelsNum + additionalPixelsYBottom + additionalPixelsYTop;
        % mean-padded photons
        photons = zeros(newYpixelsNum, newXpixelsNum, size(photons,3));
        meanPadValues = originalPhotons(1,1,:);
        xOffset = additionalPixelsXLeft;
        yOffset = additionalPixelsYBottom;
        for k = 1:size(photons,3)
            photons(:,:,k) = meanPadValues(k);
            photons(yOffset+(1:originalYpixelsNum), xOffset+(1:originalXpixelsNum),k) = originalPhotons(:,:,k);
        end
    end
    clear 'originalPhotons';
    
    % Update oiSize
    oiSize = size(photons);
    oiSize = oiSize(1:2);
    
    % Reshape the photons for efficient computations
    [photons, oiRowsNum, oiColsNum] = RGB2XWFormat(photons);
    
    % Allocate memory for noiseFreeAbsorptionsCount
    nConesNum = size(obj.coneRFpositionsMicrons,1);
    noiseFreeAbsorptionsCount = zeros(nTrials, nTimePoints, nConesNum);
        
    if (isempty(obj.fixEMobj)) || (all(emPathsMicrons(:)==0))
        % No emPath, so compute a single shot
        % Compute boost factors by which photons have to be multiplied so as
        % the account for MP density decrease with ecc
        macularPigmentDensityBoostFactors = ...
            updateMPBoostFactorsForCurrentEMpos(obj, [0 0], oiPositionsDegs, oiWave, oiSize, oiResMicrons);
    
        % Compute density of cone absosprions, by integrating photons over
        % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
        absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);

        % Compute absorptions
        %t1 = clock;
        noiseFreeAbsorptionsCount(1,1,:) = obj.integrationTime *  ...
                obj.computeAbsorptionRate(...
                    emPathsMicrons(1,1,:), ...
                    oiPositionsVectorsMicrons, ...
                    absorptionsDensityFullMap, ...
                    oiResMicrons, coneApertureDiametersMicrons, ...
                    coneIndicesInZones);
        %fprintf('Computing mean absorptions count took %f seconds.\n', etime(clock, t1));
        %t2 = clock;
        % Replicate single shot for all trials and time points
        for iTrial = 2:nTrials
            %fprintf('Replicating trial %d from first trial.\n', iTrial);
            noiseFreeAbsorptionsCount(iTrial, :, :) = noiseFreeAbsorptionsCount(1, :, :);
        end
       %fprintf('Replicating mean response for %d trials count took %f seconds.\n', nTrials, etime(clock, t2));
        
    else 
        if (~obj.eccVaryingMacularPigmentDensityDynamic)
            % No dynamic adjustment of MP density, so compute MP boost factors for the mean  eye movement position
            % across all time points and all instances. Alternatively, we
            % could do mean over all time points separately for each instance
            meanEMPosDegs = mean(mean(emPathsDegs,1),2);
            macularPigmentDensityBoostFactors = ...
                            updateMPBoostFactorsForCurrentEMpos(obj, [meanEMPosDegs(1) meanEMPosDegs(2)], oiPositionsDegs, oiWave, oiSize, oiResMicrons);
        end
        
        % Compute for emPath
        for iTrial = 1:nTrials
            if (replicateResponseToFirstEMpath) && (iTrial > 1)
                % Ideantical emPaths across trials, differing just in noise
                %fprintf('Replicating response (%d) to first trial.\n', iTrial);
                noiseFreeAbsorptionsCount(iTrial, :, :) = noiseFreeAbsorptionsCount(1, :, :);
            else
                % Different emPath for each trial
                for timePoint = 1:nTimePoints
                    
                    fprintf('\n(%s) Computing absorptions count for %d/%d time point ...', datestr(now), timePoint, nTimePoints);
                    %t1 = clock;
                    
                    if (obj.eccVaryingMacularPigmentDensityDynamic)
                        % Recompute MP boost factors for current eye movement position
                        currentEMposDegs = [emPathsDegs(iTrial, timePoint,1) emPathsDegs(iTrial, timePoint,2)];
                        macularPigmentDensityBoostFactors = ...
                            updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiWave, oiSize, oiResMicrons);
                    end
   
                    %fprintf('Computed MPBoostFactors in %2.2f seconds\n', etime(clock, t1));
                    
                    % Compute density of cone absosprions, by integrating photons over
                    % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
                    absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);
                    
                    %t1 = clock;
                    % Compute absorptions
                    noiseFreeAbsorptionsCount(iTrial, timePoint, :) = obj.integrationTime * ...
                        obj.computeAbsorptionRate(...
                        emPathsMicrons(iTrial, timePoint,:), ...
                        oiPositionsVectorsMicrons, ...
                        absorptionsDensityFullMap, ...
                        oiResMicrons, coneApertureDiametersMicrons, ...
                        coneIndicesInZones);
                    %fprintf('Computed absorptionRate in %2.2f seconds\n', etime(clock, t1));
                    
                    %fprintf(' in %f seconds.\n', etime(clock, t1));
                end % timePoint
            end
        end % iTrial
    end
    %fprintf('Tile lapsed to compute mean response: %2.2f seconds\n', etime(clock, tStart));


    responseTemporalSupport = (0:(size(noiseFreeAbsorptionsCount,2)-1)) * obj.integrationTime;
    
    if (strcmp(obj.noiseFlag, 'none'))
        noisyAbsorptionInstances = [];
    else
        % Add photon noise here.
        if (isempty(noiseSeed))
            noisyAbsorptionInstances = cMosaic.noisyInstances(noiseFreeAbsorptionsCount);
        else
            noisyAbsorptionInstances = cMosaic.noisyInstances(noiseFreeAbsorptionsCount, 'seed', noiseSeed);
        end
        %fprintf('Tile lapsed to compute Poisson noise for %d trials: %2.2f seconds\n', nTrials, etime(clock, tStart));
    end
    
    % If we have no eye movements, just return the first noise free absorptions count
    if ((isempty(obj.fixEMobj)) || (all(emPathsMicrons(:)==0)) && (nTrials > 1))
        %fprintf('Returning only the first noise free absorptions response, because there is no eye movement data\n');
        noiseFreeAbsorptionsCount = noiseFreeAbsorptionsCount(1,:,:);
    end
    
    photoCurrents = [];
    photoCurrentInstances = [];

    % All done. Restore original values
    obj.wave = originalValues.wave;
end


function macularPigmentDensityBoostFactors = updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiWave, oiSize, oiResMicrons)
    
    if (obj.eccVaryingMacularPigmentDensity)
        % Separate boost factors for each oiPixel
        macularPigmentDensityBoostFactors = obj.computeMPBoostFactors(oiPositionsDegs, currentEMposDegs, oiWave, oiSize, oiResMicrons);
    else
        % Single boost factor for the center of the mosaic
        macularPigmentDensityBoostFactor = obj.computeMPBoostFactors(obj.eccentricityDegs, currentEMposDegs, oiWave, oiSize, oiResMicrons);
        % Replicate for all oiPixels
        macularPigmentDensityBoostFactors = bsxfun(@plus, zeros(size(oiPositionsDegs,1), numel(macularPigmentDensityBoostFactor)), macularPigmentDensityBoostFactor);
    end
    
end

% Method for the validation and decoding for optional argument 'opticalImagePositionDegs'
function opticalImagePositionMicrons = validateAndDecodeOpticalImagePosition(obj, opticalImagePositionDegs)
    if ((ischar(opticalImagePositionDegs))&&(strcmp(opticalImagePositionDegs, 'mosaic-centered')))
        opticalImagePositionMicrons = obj.eccentricityMicrons;
    elseif ((isnumeric(opticalImagePositionDegs))&&(numel(opticalImagePositionDegs)==2))
        if (~isempty(obj.micronsPerDegreeApproximation))
            opticalImagePositionMicrons = opticalImagePositionDegs * obj.micronsPerDegreeApproximation;
        else
            opticalImagePositionMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(opticalImagePositionDegs);
        end
    else
        error('''opticalImagePositionDegs'' must be set to either ''mosaic-centered'' or to a 2-element vector.');
    end
end
    
% Method for the validation and decoding for optional argument 'withFixationalEyeMovements'
function [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
        validateAndDecodeFixationalEyeMovements(obj, withFixationalEyeMovements, nTrials, nTimePoints)
    
     if (withFixationalEyeMovements) && (isempty(obj.fixEMobj))
        error('cMosaic.emGenSequence() has not been called yet to generate eye movements.');
     end
    
     % emPaths/nTrials validation
    replicateResponseToFirstEMpath = false;
    if (withFixationalEyeMovements)
        % Full blown simulation. This is going to be the slowest scenario.
        % Extract the emPaths
        emPathsMicrons = obj.fixEMobj.emPosMicrons;
    	emPathsDegs = obj.fixEMobj.emPosArcMin/60;
        assert(size(emPathsMicrons,3) == 2, sprintf('The third dimension of an emPath must be 2.'));
        if (~isempty(nTrials))
           % fprintf(2,'Overiding nTrials from the attached ''fixationalEyeMovementsOBJ'' (%d). Will generate %d response instances all driven by the first emPath', ...
           %     size(emPathsMicrons,1), nTrials);
            emPathsMicrons = repmat(emPathsMicrons(1,:,:), [nTrials 1 1]);
            emPathsDegs = repmat(emPathsDegs(1,:,:), [nTrials 1 1]);
            replicateResponseToFirstEMpath = true;
        end
        if (~isempty(nTimePoints))
            error('Cannot set''nTimePoints'' when  ''withFixationalEyeMovements'' is set to true.');
        end
        nTrials = size(emPathsMicrons,1);
        nTimePoints = size(emPathsMicrons,2);
    else
        % No emPaths passed
        if (isempty(nTimePoints))
            nTimePoints = 1;
        end
        if (isempty(nTrials))
            nTrials = 1;
        end
        emPathsMicrons = zeros(nTrials,nTimePoints,2);
        emPathsDegs = zeros(nTrials,nTimePoints,2);
    end
end