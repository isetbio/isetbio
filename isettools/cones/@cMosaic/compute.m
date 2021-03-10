function [noiseFreeAbsorptionsCount, noisyAbsorptionInstances, photoCurrents, photoCurrentInstances] = compute(obj, oi, varargin)
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
    p.addParameter('nTimePoints', [], @isscalar);
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('seed', 1, @isnumeric);
    p.addParameter('verbosityLevel', 'none', @(x)ismember(x, {'default', 'min', 'max'}));
    p.parse(varargin{:});
    
    %verbosityLevel = p.Results.verbosityLevel;
    
    %tStart = clock();
    
    % Validate optional inputs
    if (p.Results.withFixationalEyeMovements) && (isempty(obj.fixEMobj))
        error('cMosaic.emGenSequence() has not been called yet to generate eye movements.');
    end
    
    
    nTimePoints = p.Results.nTimePoints;
    nTrials = p.Results.nTrials;
    noiseSeed = p.Results.seed;
    
    % emPaths/nTrials validation
    replicateResponseToFirstEMpath = false;
    if (p.Results.withFixationalEyeMovements)
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
    
    % Retrieve oiRes and oiSize
    oiSize  = oiGet(oi, 'size');
    oiResMicrons = oiGet(oi, 'height spatial resolution')*1e6;
    
    % Generate oiPositions
    oiYPosMicrons = (1:oiSize(1))*oiResMicrons;
    oiXPosMicrons = (1:oiSize(2))*oiResMicrons;
    oiXPosMicrons = oiXPosMicrons - mean(oiXPosMicrons);
    oiYPosMicrons = oiYPosMicrons - mean(oiYPosMicrons);
    
    minEMpos = squeeze(min(emPathsMicrons,[],1));
    maxEMpos = squeeze(max(emPathsMicrons,[],1));
    
    if (obj.minRFpositionMicrons(1)+minEMpos(1) < min(oiXPosMicrons))
        fprintf(2,'Left side of mosaic extends beyond the optical image. \nExpect artifacts there. Increase optical image size to avoid these.\n');
    end
    if (obj.minRFpositionMicrons(2)+minEMpos(2) < min(oiYPosMicrons))
        fprintf(2,'Bottom side of mosaic extends beyond the optical image. \nExpect artifacts there. Increase optical image size to avoid these.\n');
    end
    
    
    if (obj.maxRFpositionMicrons(1)+maxEMpos(1) > max(oiXPosMicrons))
        fprintf(2,'Right side of mosaic extends beyond the optical image. \nExpect artifacts there. Increase optical image size to avoid these.\n');
    end
    if (obj.maxRFpositionMicrons(2)+maxEMpos(2) > max(oiYPosMicrons))
        fprintf(2, 'Top side of mosaic extends beyond the optical image. \nExpect artifacts there. Increase optical image size to avoid these.\n');
    end

    
    if (~isempty(obj.micronsPerDegreeApproximation))
        oiXPosDegrees = oiXPosMicrons/obj.micronsPerDegreeApproximation;  
        oiYPosDegrees = oiYPosMicrons/obj.micronsPerDegreeApproximation;
    else
        oiXPosDegrees = RGCmodels.Watson.convert.rhoMMsToDegs(oiXPosMicrons*1e-3);
        oiYPosDegrees = RGCmodels.Watson.convert.rhoMMsToDegs(oiYPosMicrons*1e-3);
    end

    [oiPositionsDegsX, oiPositionsDegsY] = meshgrid(oiXPosDegrees, oiYPosDegrees);
    oiPositionsDegs = [oiPositionsDegsX(:), oiPositionsDegsY(:)];
    
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
   
    % Reshape the photons for efficient computations
    [photons, oiRowsNum, oiColsNum] = RGB2XWFormat(photons);
    
    % Compute boost factors by which photons have to be multiplied so as
    % the account for MP density decrease with ecc
    %t1 = clock;
    currentEMposDegs = [emPathsDegs(1, 1,1) emPathsDegs(1,1,2)];
    macularPigmentDensityBoostFactors = ...
        updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiSize, oiResMicrons);
    
    %fprintf('Computing ecc-based MP boosting factors took %f seconds.\n', etime(clock, t1));
        
    % Allocate memory for noiseFreeAbsorptionsCount
    nConesNum = size(obj.coneRFpositionsMicrons,1);
    noiseFreeAbsorptionsCount = zeros(nTrials, nTimePoints, nConesNum);
        
    if (isempty(obj.fixEMobj)) || (all(emPathsMicrons(:)==0))
        % No emPath, so compute a single shot
        
        % Compute density of cone absosprions, by integrating photons over
        % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
        absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);

        % Compute absorptions
        %t1 = clock;
        noiseFreeAbsorptionsCount(1,1,:) = obj.integrationTime *  ...
                obj.computeAbsorptionRate(...
                    emPathsMicrons(1,1,:), ...
                    [oiXPosMicrons(:) oiYPosMicrons(:)], ...
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
                    
                    % Compute boost factors by which photons have to be multiplied so as
                    % the account for MP density decrease with ecc
                    if (obj.eccVaryingMacularPigmentDensityDynamic)
                        currentEMposDegs = [emPathsDegs(iTrial, timePoint,1) emPathsDegs(iTrial, timePoint,2)];
                        macularPigmentDensityBoostFactors = ...
                            updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiSize, oiResMicrons);
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
                        [oiXPosMicrons(:) oiYPosMicrons(:)], ...
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


function macularPigmentDensityBoostFactors = updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiSize, oiResMicrons)
    if (obj.eccVaryingMacularPigmentDensity)
        % Separate boost factors for each oiPixel
        macularPigmentDensityBoostFactors = obj.computeMPBoostFactors(oiPositionsDegs, currentEMposDegs, oiSize, oiResMicrons);
    else
        % Single boost factor for the center of the mosaic
        macularPigmentDensityBoostFactor = obj.computeMPBoostFactors(obj.eccentricityDegs, currentEMposDegs, oiSize, oiResMicrons);
        % Replicate for all oiPixels
        macularPigmentDensityBoostFactors = bsxfun(@plus, zeros(size(oiPositionsDegs,1), numel(macularPigmentDensityBoostFactor)), macularPigmentDensityBoostFactor);
    end
    
end

