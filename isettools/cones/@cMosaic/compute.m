function [noiseFreeAbsorptionsCount, noisyAbsorptionInstances, ...
    photoCurrents, photoCurrentInstances, responseTemporalSupport] = ...
    compute(obj, oi, varargin)
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
%    noiseFreeAbsorptionsCount - [nTrials, nTimePoints, nConesNum] matrix
%                                 with noise-free cone excitations
%    noisyAbsorptionInstances - [nTrials, nTimePoints, nConesNum] matrix
%                                 with noisy cone excitation instances

% Optional key/value pairs:
%    'emPaths'           - Eye movement paths. Either empty or a matrix of [nTrials x N x 2]. Default is [].
%                          

    % Parse input
    p = inputParser;
    p.addRequired('oi', @(x)((isa(x, 'struct')) || (isa(x, 'oiSequence')) || (isa(x, 'oiArbitrarySequence'))))
    p.addParameter('withFixationalEyeMovements', false, @islogical);
    p.addParameter('opticalImagePositionDegs', obj.opticalImagePositionDegs, @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
    p.addParameter('nTimePoints', [], @isscalar);
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('seed', 1, @isnumeric);
    p.addParameter('verbosityLevel', 'none', @(x)ismember(x, {'default', 'min', 'max'}));
    p.parse(oi,varargin{:});
    

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
    if (isa(oi, 'struct'))
        [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
            validateAndDecodeFixationalEyeMovementsForSingleOI(obj, p.Results.withFixationalEyeMovements, nTrials, nTimePoints);
    else
        [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
            validateAndDecodeFixationalEyeMovementsForOISequence(obj, p.Results.withFixationalEyeMovements, nTrials, nTimePoints, oi);
    end
    
    % Decode the opticalImagePositionDegs optional input
    opticalImagePositionMicrons = validateAndDecodeOpticalImagePosition(obj,p.Results.opticalImagePositionDegs);
    
    % Match the wavelength support of the object to that of the input optical image
    % First save the current value so we can restore it once we are done
    % with this computation
    originalValues.wave = obj.wave;
    
    if (isa(oi, 'struct'))
        % Set the new wavelength. The cMosaic object has listeners on wave
        % that set the wavelength support of the attached pigment and macular
        % properties
        obj.wave = oiGet(oi, 'wave');

        % Compute wavelength-spacing scaled quantal efficiencies
        scaledQE = obj.qe * oiGet(oi, 'bin width');

        % Retrieve oiRes
        oiResMicrons = oiGet(oi, 'height spatial resolution')*1e6;

        % Retrieve wavelength support
        oiWave = oiGet(oi, 'wave');
        
        % Retrieve spatial support in meters
        spatialSupportMeters = oiGet(oi, 'spatial support');
    else
        % An oiSequence, so do as above using the first OI frame
        firstOI = oi.frameAtIndex(1);
        obj.wave = oiGet(firstOI, 'wave');
        scaledQE = obj.qe * oiGet(firstOI, 'bin width');
        oiResMicrons = oiGet(firstOI, 'height spatial resolution')*1e6;
        oiWave = oiGet(firstOI, 'wave');
        spatialSupportMeters = oiGet(firstOI, 'spatial support');
    end
    
    % Generate oiPositions
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
        additionalPixels.xLeft = round((min(spatialSupportXMicrons) - (obj.minRFpositionMicrons(1)+minEMpos(1)))/dx);
    else
        additionalPixels.xLeft = 0;
    end
    
    if (obj.minRFpositionMicrons(2)+minEMpos(2) < min(spatialSupportYMicrons))
        additionalPixels.yBottom = round((min(spatialSupportYMicrons) - (obj.minRFpositionMicrons(2)+minEMpos(2)))/dy);
    else
        additionalPixels.yBottom = 0;
    end
    
    if (obj.maxRFpositionMicrons(1)+maxEMpos(1) > max(spatialSupportXMicrons))
        additionalPixels.xRight = round((obj.maxRFpositionMicrons(1)+maxEMpos(1) - max(spatialSupportXMicrons))/dx);
    else
        additionalPixels.xRight = 0;
    end
    
    if (obj.maxRFpositionMicrons(2)+maxEMpos(2) > max(spatialSupportYMicrons))
        additionalPixels.yTop = round((obj.maxRFpositionMicrons(2)+maxEMpos(2) - max(spatialSupportYMicrons))/dy);
    else
        additionalPixels.yTop = 0;
    end
    
    % Extend spatial support vectors as needed
    if (additionalPixels.xLeft > 0)
        extraXPos = -fliplr((1:additionalPixels.xLeft)*dx);
        spatialSupportXMicrons = cat(2, extraXPos+spatialSupportXMicrons(1), spatialSupportXMicrons); 
    end
    if (additionalPixels.xRight > 0)
        extraXPos = (1:additionalPixels.xRight)*dx;
        spatialSupportXMicrons = cat(2, spatialSupportXMicrons, extraXPos+spatialSupportXMicrons(end)); 
    end
    
    if (additionalPixels.yBottom > 0)
        extraYPos = -fliplr((1:additionalPixels.yBottom)*dy);
        spatialSupportYMicrons = cat(1, extraYPos'+spatialSupportYMicrons(1), spatialSupportYMicrons);
    end
    if (additionalPixels.yTop > 0)
        extraYPos = (1:additionalPixels.yTop)*dy;
        spatialSupportYMicrons = cat(1,  spatialSupportYMicrons, extraYPos'+spatialSupportYMicrons(end));
    end
    
    spatialSupportXDegrees = obj.distanceMicronsToDistanceDegreesForCmosaic(spatialSupportXMicrons);
    spatialSupportYDegrees = obj.distanceMicronsToDistanceDegreesForCmosaic(spatialSupportYMicrons);
    
    
    [oiPositionsDegsXgrid, oiPositionsDegsYgrid] = meshgrid(spatialSupportXDegrees, spatialSupportYDegrees);
    oiPositionsDegs = [oiPositionsDegsXgrid(:), oiPositionsDegsYgrid(:)];
    oiPositionsVectorsMicrons = {spatialSupportYMicrons(:), spatialSupportXMicrons(:)};
    

    % Allocate memory for noiseFreeAbsorptionsCount
    nConesNum = size(obj.coneRFpositionsMicrons,1);
    noiseFreeAbsorptionsCount = zeros(nTrials, nTimePoints, nConesNum);
    
    % Form responseTemporalSupport
    responseTemporalSupport = (0:(size(noiseFreeAbsorptionsCount,2)-1)) * obj.integrationTime;
    
    % Parse time
    if (isa(oi, 'struct'))
        oiTimeAxis = [0];
    else
        oiTimeAxis = oi.timeAxis;
    end
    oiFramesNum = numel(oiTimeAxis);
    
            
    for oiFrame = 1:oiFramesNum
        
        if (isa(oi, 'struct'))
            % Retrieve retinal irradiance in photons
            photons = oiGet(oi, 'photons');

        else
            % Retrieve retinal irradiance in photons
            photons = oiGet(oi.frameAtIndex(oiFrame), 'photons');
        end


        % Flip the optical image upside-down because the y-coords in the
        % oi spatial support vectors increase from top -> bottom (y-coords in an image)
        % whereas cone y-positions increase from bottom -> top 
        for k = 1:size(photons,3)
            photons(:,:,k) = flipud(squeeze(photons(:,:,k)));
        end

        % Phad photons field as needed to accomodate eye movements
        photons = padPhotons(photons, additionalPixels);
    
        % Update oiSize
        oiSize = size(photons);
        oiSize = oiSize(1:2);

        % Reshape the photons for efficient computations
        [photons, oiRowsNum, oiColsNum] = RGB2XWFormat(photons);

        % Single OI
        if (isa(oi, 'struct'))
            
            if (isempty(obj.fixEMobj)) || (all(emPathsMicrons(:)==0))
                % No emPath, so compute a single shot
                %
                % Compute boost factors by which photons have to be
                % multiplied so as to account for MP density decrease with
                % eccentricity
                macularPigmentDensityBoostFactors = ...
                    updateMPBoostFactorsForCurrentEMpos(obj, [0 0], oiPositionsDegs, oiWave, oiSize, oiResMicrons);
                
                % Compute density of cone absosprions, by integrating photons over
                % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
                absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);
                 
                % Compute absorptions
                noiseFreeAbsorptionsCount(1,1,:) = obj.integrationTime *  ...
                    obj.computeAbsorptionRate(emPathsMicrons(1,1,:), ...
                                              oiPositionsVectorsMicrons, ...
                                              absorptionsDensityFullMap, ...
                                              oiResMicrons);

                % Replicate single shot for all trials and time points
                for timePoint = 1:nTimePoints
                    noiseFreeAbsorptionsCount(1, timePoint, :) = noiseFreeAbsorptionsCount(1,1,:);
                end
                
                for iTrial = 2:nTrials
                    %fprintf('Replicating trial %d from first trial.\n', iTrial);
                    noiseFreeAbsorptionsCount(iTrial, :, :) = noiseFreeAbsorptionsCount(1, :, :);
                end
                
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
                            
                            if (obj.eccVaryingMacularPigmentDensityDynamic)
                                % Recompute MP boost factors for current eye movement position
                                currentEMposDegs = [emPathsDegs(iTrial, timePoint,1) emPathsDegs(iTrial, timePoint,2)];
                                macularPigmentDensityBoostFactors = ...
                                    updateMPBoostFactorsForCurrentEMpos(obj, currentEMposDegs, oiPositionsDegs, oiWave, oiSize, oiResMicrons);
                            end
                            
                            % Compute density of cone absosprions, by integrating photons over
                            % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
                            absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);
                            
                            % Compute absorptions
                            noiseFreeAbsorptionsCount(iTrial, timePoint, :) = obj.integrationTime * ...
                                obj.computeAbsorptionRate(emPathsMicrons(iTrial, timePoint,:), ...
                                                          oiPositionsVectorsMicrons, ...
                                                          absorptionsDensityFullMap, ...
                                                          oiResMicrons);
                            
                        end % timePoint
                    end
                end % iTrial
                
            end
            
        % OISequence
        else
            if (isempty(obj.fixEMobj)) || (all(emPathsMicrons(:)==0))
                % No emPath
                
                if (oiFrame == 1)
                    % Compute boost factors by which photons have to be multiplied so as
                    % the account for MP density decrease with ecc. This
                    % needs to be done only once (we do it for the first frame only)
                    macularPigmentDensityBoostFactors = ...
                        updateMPBoostFactorsForCurrentEMpos(obj, [0 0], oiPositionsDegs, oiWave, oiSize, oiResMicrons);
                end
                
                % Compute density of cone absosprions, by integrating photons over
                % wavelength. The size of abosrptionsDensity is [oiRows x oiCols x coneTypes]
                absorptionsDensityFullMap = XW2RGBFormat((photons .* macularPigmentDensityBoostFactors) * scaledQE, oiRowsNum, oiColsNum);  
                
                % Compute absorptions for this frame
                noiseFreeAbsorptionsCount(1,oiFrame,:) = obj.integrationTime *  ...
                    obj.computeAbsorptionRate(emPathsMicrons(1,1,:), ...
                                              oiPositionsVectorsMicrons, ...
                                              absorptionsDensityFullMap, ...
                                              oiResMicrons);

                % Replicate this frame for all trials
                for iTrial = 2:nTrials
                    noiseFreeAbsorptionsCount(iTrial, oiFrame, :) = noiseFreeAbsorptionsCount(1,oiFrame,:);
                end
            end
            
        end  % oiSequence
    end % for oiFrame
    
    
    % Perform cone coupling at the level of cone excitations, if specified
    if (~isempty(obj.coneCouplingLambda)) && (abs(obj.coneCouplingLambda)>0)
        fprintf(2, 'Electrically coupling cone responses at the level of CONE EXCITATIONS with lambda factor = %2.3f\n', obj.coneCouplingLambda);
        noiseFreeAbsorptionsCount = obj.electricallyCoupleConeResponses(noiseFreeAbsorptionsCount);
    end
    
    
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

    % Save last absorptionsDensityFullMap, oiPositionsVectorsMicrons
    obj.absorptionsDensityFullMap = absorptionsDensityFullMap;
    obj.absorptionsDensitySpatialSupportMicrons = oiPositionsVectorsMicrons;
    
    % All done. Restore original values
    obj.wave = originalValues.wave;
end


function photons = padPhotons(photons, additionalPixels)
    if ((additionalPixels.xLeft > 0) || (additionalPixels.xRight > 0) || ...
        (additionalPixels.yBottom > 0) || (additionalPixels.yTop > 0))
        originalPhotons = photons;
        originalXpixelsNum = size(photons,2);
        originalYpixelsNum = size(photons,1);
        newXpixelsNum = originalXpixelsNum + additionalPixels.xLeft + additionalPixels.xRight;
        newYpixelsNum = originalYpixelsNum + additionalPixels.yBottom + additionalPixels.yTop;
        % mean-padded photons
        photons = zeros(newYpixelsNum, newXpixelsNum, size(photons,3));
        meanPadValues = originalPhotons(1,1,:);
        xOffset = additionalPixels.xLeft;
        yOffset = additionalPixels.yBottom;
        for k = 1:size(photons,3)
            photons(:,:,k) = meanPadValues(k);
            photons(yOffset+(1:originalYpixelsNum), xOffset+(1:originalXpixelsNum),k) = originalPhotons(:,:,k);
        end
    end
    clear 'originalPhotons';
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
        opticalImagePositionMicrons = obj.distanceDegreesToDistanceMicronsForCmosaic(opticalImagePositionDegs);
    else
        error('''opticalImagePositionDegs'' must be set to either ''mosaic-centered'' or to a 2-element vector.');
    end
end
    

% Method for the validation and decoding for optional argument
% 'withFixationalEyeMovements' when we are dealing with an OIsequence
function [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
            validateAndDecodeFixationalEyeMovementsForOISequence(obj, withFixationalEyeMovements, nTrials, nTimePoints, oiSequence)
        
    if (withFixationalEyeMovements) && (isempty(obj.fixEMobj))
        error('cMosaic.emGenSequence() has not been called yet to generate eye movements.');
    end
     
    % Get the time axis of the oiSequence
    timeAxis = oiSequence.timeAxis;
        
    replicateResponseToFirstEMpath = false;
    if (withFixationalEyeMovements)
        error('cMosaic.compute(): We have not yet implemented this method for fixationalEyeMovements with an OIsequence');
    else
        if (~isempty(nTimePoints))
            fprintf('cMosaic.compute() with oiSequence: ignoring ''nTimePoints'' (set to %d) parameter, and setting it to the length of the oiSequence.', nTimePoints);
        end
        % No emPaths passed: make the nTimePoints = legth of oiSequence
        nTimePoints = numel(timeAxis);
        % Also ensure that the cone mosaic's integrationTime matches the
        % oiSequence frame duration
        if (numel(timeAxis) > 1)
            assert(obj.integrationTime == timeAxis(2)-timeAxis(1), ...
                'cMosaic.compute() with oiSequence: the mosaic''s integrationTime (%2.2f msec) does not match the oiSequence.frame duration (%2.2f msec).', ...
                obj.integrationTime*1000, (timeAxis(2)-timeAxis(1))*1000)
        end
        
        if (isempty(nTrials))
            nTrials = 1;
        end
        emPathsMicrons = zeros(nTrials,nTimePoints,2);
        emPathsDegs = zeros(nTrials,nTimePoints,2);
    end
    
end


        
% Method for the validation and decoding for optional argument
% 'withFixationalEyeMovements' when we are dealing with a singleOI
function [emPathsDegs, emPathsMicrons, nTrials, nTimePoints, replicateResponseToFirstEMpath] = ...
        validateAndDecodeFixationalEyeMovementsForSingleOI(obj, withFixationalEyeMovements, nTrials, nTimePoints)
    
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