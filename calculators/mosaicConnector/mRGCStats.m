function [mRGCSpacingInMicrons, conesToRGCratio] = mRGCStats(rfPositions, nSamples, whichEye)

    % Find range of retinal positions in microns that we need to compute density for
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
    conesToMRGCratio2Dmap(conesToMRGCratio2Dmap<1) = 1;
    
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