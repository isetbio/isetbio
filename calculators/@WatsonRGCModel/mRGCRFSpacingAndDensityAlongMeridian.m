function [mRGCRFSpacing, mRGCRFDensity, rightEyeRetinalMeridianName] = ...
    mRGCRFSpacingAndDensityAlongMeridian(obj, eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits, varargin)

% Input
%   eccentricities      1-D vector with eccentricities (specified in eccUnits)
%   rightEyeVisualFieldMeridianName        name of the meridian in Watson's reference (visual
%                       field of the right eye, with temporal being at 0 degs
%                       superior at 90, nasal at 18, & inferior at 270 degs

    % Parse input
    p = inputParser;       
    p.addParameter('adjustForISETBioConeDensity', false, @islogical);
    p.addParameter('subtype', 'ONOFF', @(x)(ismember(x, {'ON', 'OFF', 'ONOFF'})));
    p.parse(varargin{:});
   
    % Validate eccUnits
    obj.validateEccUnits(eccUnits);
    
    % Validate densityUnits
    obj.validateDensityUnits(densityUnits);
    
    % Make sure eccentricities is a 1xN vector
    if (size(eccentricities,1)>1)
        eccentricities = eccentricities';
    end
    assert(size(eccentricities,1) == 1, 'Eccentricities must be a 1xN vector');
    
    % Convert passed eccentricities as requested
    switch (eccUnits)
        case obj.visualDegsEccUnits
            % Convert ecc from degs to retinal MMs
            eccMM = obj.rhoDegsToMMs(eccentricities);
            eccDegs = eccentricities;
        case obj.retinalMMEccUnits
            eccMM = eccentricities;
            eccDegs = obj.rhoMMsToDegs(eccMM);
    end
    
    
    % Compute total RGC RF density along the requested meridian
    totalRGCRFdensityAlongMeridian = obj.totalRGCRFDensityAlongMeridian(...
        eccDegs, rightEyeVisualFieldMeridianName, obj.visualDegsEccUnits, densityUnits);
    
    % Retrieve percentage of total ganglion cells that are midget at 0 eccentricity
    percentageOfMidgetRGCsAtZeroEccentricity = obj.f0;
    
    % mRGC fraction (ot total RGCs) along meridian. This is equation (7) in the Watson (2014) paper.
    midgetRGCFractionAlongMeridian = percentageOfMidgetRGCsAtZeroEccentricity ./ (1 + eccDegs/41.03);
   
    % Compute mRGC RF density along the requested meridian. This is equation (8) in the Watson (2014) paper.
    mRGCRFDensity = midgetRGCFractionAlongMeridian .* totalRGCRFdensityAlongMeridian;
    
    switch(p.Results.subtype)
        case 'ON'
            mRGCRFDensity = 0.5*mRGCRFDensity;
        case 'OFF'
            mRGCRFDensity = 0.5*mRGCRFDensity;
        case 'ONOFF'
            % do nothing
        otherwise
            error('subtype must be ''ON'', ''OFF'', or ''ONOFF''.');
    end
    
    
    if (p.Results.adjustForISETBioConeDensity)
        [~, coneDensityWatson] = obj.coneRFSpacingAndDensityAlongMeridian( ...
        	eccDegs, rightEyeVisualFieldMeridianName, obj.visualDegsEccUnits, densityUnits, ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', true);
        
        [~, coneDensityISETBio] = obj.coneRFSpacingAndDensityAlongMeridian( ...
        	eccDegs, rightEyeVisualFieldMeridianName, obj.visualDegsEccUnits, densityUnits, ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    
        mRGCToConeRatio = mRGCRFDensity ./ coneDensityWatson;
        mRGCRFDensity = coneDensityISETBio .* mRGCToConeRatio;
    end
    
    % Compute mRGC RF spacing from their density. 
    mRGCRFSpacing = obj.spacingFromDensity(mRGCRFDensity);
    
    % Return the name of the corresponding retinal meridian.
    [~, ~, rightEyeRetinalMeridianName] = obj.isetbioRetinalAngleForWatsonMeridian(rightEyeVisualFieldMeridianName);
    
end


