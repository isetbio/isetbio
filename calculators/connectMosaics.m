function connectMosaics()
   
    micronsPerDeg = 300;
    whichEye = 'right';
    mosaicFOVDegs = 30;
    eccentricitySamplesNum = 48;

    micronsPerDeg = 300;
    analyzedRadiusDeg = 2.5;
    workingRadiusMicrons = analyzedRadiusDeg*micronsPerDeg;

    
    neuronalType = 'mRGC';
    RGCRFPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum);
    [RGCRFSpacings, conesToRGCratios] = mRGCStats(RGCRFPositions, 128, whichEye);

    neuronalType = 'cone';
    conePositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum);
    
    %Only keep cones and RGC within the working radius
    eccCones = sqrt(sum(conePositions.^2,2));
    idx = find(eccCones<workingRadiusMicrons);
    conePositions = conePositions(idx,:);
    
    eccRGC = sqrt(sum(RGCRFPositions.^2,2));
    idx = find(eccRGC<workingRadiusMicrons);
    RGCRFPositions = RGCRFPositions(idx,:);
    RGCRFSpacings = RGCRFSpacings(idx);
    conesToRGCratios = conesToRGCratios(idx);

    
    fprintf('Loaded %2.0f RGCs\n', size(RGCRFPositions,1));
    fprintf('Loaded %2.0f cones\n', size(conePositions,1));
    

    connectionMatrix = matchMosaics(conePositions, RGCRFPositions, RGCRFSpacings, conesToRGCratios);

    maxConeWeights = max(connectionMatrix,[],2);
    connectionMatrix = bsxfun(@times, connectionMatrix, 1./maxConeWeights);
    
    roiCenter = [1 1]*300;
    roiSize = [0.3 0.3]*300;
    mRGCindices = find(...
        (abs(RGCRFPositions(:,1)-roiCenter(1))< roiSize(1)/2) & ...
        (abs(RGCRFPositions(:,2)-roiCenter(2))< roiSize(2)/2) );
    coneIndices = find(...
        (abs(conePositions(:,1)-roiCenter(1))< roiSize(1)/2) & ...
        (abs(conePositions(:,2)-roiCenter(2))< roiSize(2)/2) );
    
    figure(1); clf;
    hold on;
    % Plot connections
    for m = 1:numel(mRGCindices)
        mRGCindex = mRGCindices(m);
        
        mPos = RGCRFPositions(mRGCindex,:);
        
        plot(mPos(1)+cosd(0:30:360)*RGCRFSpacings(mRGCindex)/2, ...
             mPos(2)+sind(0:30:360)*RGCRFSpacings(mRGCindex)/2, 'k--', 'LineWidth', 1);
         
        for c = 1:numel(coneIndices)
            cIndex = coneIndices(c);
            cPos = conePositions(cIndex,:);
            connectionStrength = connectionMatrix(mRGCindex,cIndex);
            if (connectionStrength>0)
                plot([mPos(1) cPos(1)],[mPos(2) cPos(2)], '-', ...
                    'Color', connectionStrength*[1 0 0], 'LineWidth', connectionStrength*3);
            end
            
            drawnow
        end
    end
    %plot(RGCRFPositions(:,1), RGCRFPositions(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 14); hold on;
    plot(conePositions(coneIndices,1), conePositions(coneIndices,2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
    
    axis 'square'
    
end


function connectionMatrix = matchMosaics(conePositions, RGCRFPositions, RGCRFSpacings, conesToRGCratios)


    conesNum = size(conePositions,1);
    RGCsNum = size(RGCRFPositions,1);
    
    connectionMatrix = zeros(RGCsNum, conesNum, 'single');
    
    RGCinputsNum = zeros(RGCsNum,1);
    
    coneEcc = sqrt(sum(conePositions.^2,2));
    [~,sortedConeIndices] = sort(coneEcc, 'ascend');
    
        
    % Connect every cone to its nearest RGC 
    for cIndex = 1:conesNum
        theTargetConeIndex = sortedConeIndices(cIndex);
        theTargetConePos = conePositions(theTargetConeIndex,:);
        
        % distance to all RGCs
        distanceToAllRGCs = sqrt(sum((bsxfun(@minus, RGCRFPositions, theTargetConePos)).^2,2));
        
        % find the closest RGC 
        [~,closestRGCindex] = min(distanceToAllRGCs);
        
        % Choose the spacing of that RGC as the search radius
        searchRadius = RGCRFSpacings(closestRGCindex);
        
        % Find nearby RGCs (within a searchRadius from the target cone)
        possibleRGCindicesForConnection = find(distanceToAllRGCs <= searchRadius);
          
        % Examine which of those nearby RGCs to connect to
        if (~isempty(possibleRGCindicesForConnection))

            % Likelihood of connecting is inversely proportional to distance
            willAcceptAnotherConnection = 1./distanceToAllRGCs(possibleRGCindicesForConnection);
            willAcceptAnotherConnection(willAcceptAnotherConnection>100) = 100;
            
            % Do not connect connections to RGCs that already have a connection 
            for k = 1:numel(possibleRGCindicesForConnection)
                theRGCindex = possibleRGCindicesForConnection(k);
                if (RGCinputsNum(theRGCindex) > 0)
                    willAcceptAnotherConnection(k) = 0;
                end
            end
            
            % Find RGC with P of connecting to
            [maxAcceptance,idx] = max(willAcceptAnotherConnection);
            
            if (maxAcceptance > 0)
                matchedRGCindex = possibleRGCindicesForConnection(idx);
                %fprintf('Connecting 1st cone to RGC %d\n',  matchedRGCindex);
                connectionStrength = 1;
            else
                % All possible RGCs have at least 1 input, so find the RGC
                % with the minimal cone connections that is also closest to the target cone
                % The min number of cone inputs in this collection of RGCs
                minConeInputs = min(RGCinputsNum(possibleRGCindicesForConnection));
                
                % The RGCs which have this minimum number of inputs
                idx = find(RGCinputsNum(possibleRGCindicesForConnection) == minConeInputs);
                RGCindicesWithLeastConeInputs = possibleRGCindicesForConnection(idx);
                
                % Sort them according to their distance to the target cone
                [~,idx] = min(distanceToAllRGCs(RGCindicesWithLeastConeInputs));
                
                % Select the RGC with the closest distance
                matchedRGCindex = RGCindicesWithLeastConeInputs(idx);
                
                %fprintf('\t Connecting %dnd cone to RGC %d\n', RGCinputsNum(matchedRGCindex)+1, matchedRGCindex);
                connectionStrength = 1/minConeInputs;
            end
            
            connectionMatrix(matchedRGCindex,theTargetConeIndex) = single(connectionStrength);
                 
            % Change this RGC position to match that of the cone
            %RGCRFPositionsNew(matchedRGCindex,:) = theTargetConePos;

            % Update the input connections to this RGC
            RGCinputsNum(matchedRGCindex) = RGCinputsNum(matchedRGCindex)+1;
        end
        
    end % cIndex

end


function [RGCRFPositionsNew, connectionStrengths] = matchMosaicsDendrites(conePositions, RGCRFPositions, RGCRFSpacings)

    rgcSigma = 0.3*RGCRFSpacings/2;
    coneSigma = 4;
    
    minPos = min(conePositions,[],1);
    maxPos = max(conePositions,[],1);
    
    deltaXY = 1;
    X = minPos(1):deltaXY:maxPos(1);
    Y = minPos(2):deltaXY:maxPos(2);
    X
    Y
    pause
    [X,Y] = meshgrid(X,Y);
    
    conesNum = size(conePositions,1);
    RGCsNum = size(RGCRFPositions,1);
    
    
    for cIndex = 1:conesNum
        theTargetConePos = conePositions(cIndex,:);
        coneSpread(cIndex,:,:) = exp(-0.5*((X-theTargetConePos(1))/coneSigma).^2) .* ...
                                 exp(-0.5*((Y-theTargetConePos(2))/coneSigma).^2);
                             figure(1)
                             imagesc(squeeze(coneSpread(cIndex,:,:)))
                             drawnow;
    end
    
    for mIndex = 1:RGCsNum
        theTargetRGCPos = RGCRFPositions(mIndex,:);
        rgcSpread(mIndex,:,:) = exp(-0.5*((X-theTargetRGCPos(1))/rgcSigma(mIndex)).^2) .* ...
                                exp(-0.5*((Y-theTargetRGCPos(2))/rgcSigma(mIndex)).^2);
                            figure(1)
                             imagesc(squeeze(rgcSpread(mIndex,:,:)))
                             drawnow;
    end
    
    for mRGCindex = 1:size(RGCRFPositions,1)
        rgcDendriteRF = squeeze(rgcSpread(mIndex,:,:));
        Z = sum(sum(rgcDendriteRF.*rgcDendriteRF));
        for cIndex = 1:conesNum
            coneProcessesRF = squeeze(coneSpread(cIndex,:,:));
            tmp = sum(sum(rgcDendriteRF.*coneProcessesRF))/Z
            connectionStrengths(mRGCindex,cIndex) = tmp;
        end
        connectionStrengths(mRGCindex,:) = connectionStrengths(mRGCindex,:)/max(connectionStrengths(mRGCindex,:));
        RGCRFPositionsNew(mRGCindex,:) = RGCRFPositions(mRGCindex,:);
        
    end
    
    connectionStrengths(connectionStrengths<0.01) = 0;
end





function rfPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum)
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    
    
    RGCMosaicFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_prctile%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, 99));

    load(RGCMosaicFileName, 'rfPositionsHistory');
    rfPositions = double(squeeze(rfPositionsHistory(end,:,:)));   
end


function [mRGCSpacingInMicrons, conesToRGCratio, eccentricitiesInMicrons] = mRGCStats(rfPositions, nSamples, whichEye)

    % Find range of retinal positions in microns that we need to compute
    % density for
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    idx = find(eccentricitiesInMicrons<2.1);
    t = rfPositions(idx,:);
    tt = pdist2(t,t);
    tt(tt==0) = Inf;
    minSeparationMicrons = min(tt(:));
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Support of retinal positions in mm
    xPosMM = [0 logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples)]/1e3;
    
    switch whichEye
        case 'left'
            theView = 'left eye retina';
        case 'right'
            theView = 'right eye retina';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end

    WatsonRGCModelObj = WatsonRGCModel();
    
    % Convert retinal mm to visual degs
    eccDegs = WatsonRGCModelObj.rhoMMsToDegs(xPosMM);

    % Compute mRGC density map
    [mRGCDensity2DMap, mRGCMeridianDensities, densitySupportMM, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DmRGCRFDensity(eccDegs, theView);

    % Density is for both types of mRGCs (ON + OFF), so we need density for
    % one type, which is half (assuming equal numerosities of ON and OFF cells)
    mRGCDensity2DMap = 0.5*mRGCDensity2DMap;
    
    % Compute cone density map
    [coneDensity2DMap, coneMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DConeRFDensity(eccDegs, theView);
   
    
    % Make sure results are returned in retinal mm units
    assert((strcmp(densityUnits, WatsonRGCModelObj.retinalMMDensityUnits)) && ...
           (strcmp(supportUnits, WatsonRGCModelObj.retinalMMEccUnits)), ...
           sprintf('Expected mm units, but got ''%s'' and ''%s'' instead.', supportUnits, densityUnits)); 

    % Compute cone to mRGC ratio map
    conesToMRGCratio2Dmap = coneDensity2DMap./mRGCDensity2DMap;
    
    % Convert the density map into a spacing map
    mRGCSpacing2DMapMM = WatsonRGCModelObj.spacingFromDensity(mRGCDensity2DMap);
    
    % Convert to microns from mm
    mRGCSpacing2DMapMicrons = mRGCSpacing2DMapMM*1e3;
    densitySupportMicrons = densitySupportMM*1e3;

    % Create a scatterred interpolant function  so we can
    % compute mRGC RF spacing at the actual mRGC positions
    [X,Y] = meshgrid(squeeze(densitySupportMicrons(1,:)), squeeze(densitySupportMicrons(2,:)));
    Fspacing = scatteredInterpolant(X(:),Y(:),mRGCSpacing2DMapMicrons(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    mRGCSpacingInMicrons = Fspacing(rfPositions(:,1), rfPositions(:,2));
    
    % Create a scatterred interpolant function  so we can
    % compute cone-to-mRGC ratios at the actual mRGC positions
    Fratio = scatteredInterpolant(X(:),Y(:),conesToMRGCratio2Dmap(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    conesToRGCratio = Fratio(rfPositions(:,1), rfPositions(:,2));
    
end



