function [mRGCSpacings, mRGCRFDensities] = mRGCRFSpacingAndDensityAtRetinalPositions(obj, rfPositions, whichEye, posUnits, densityUnits, varargin)
% Return mRGCRF spacing (in mm) and densities (in units per retinal mm^2, at each of the input
% retinal positions (specified in retinal mm, for either the left or the
% right eye.
    
     % Parse input
    p = inputParser;       
    p.addParameter('adjustForISETBioConeDensity', false, @islogical);
    p.addParameter('subtype', 'ONOFF', @(x)(ismember(x, {'ON', 'OFF', 'ONOFF'})));
    p.parse(varargin{:});
    
    % Compute variation along each of the enumerated meridians
    meridianRFDensities = zeros(numel(obj.enumeratedMeridianNames),size(rfPositions,1));
   
    % Convert meridians from RE visual field to left/right retina,
    % dependin on whichEye
    mappedMeridianName = containers.Map();
    switch whichEye
        case 'left'
            %'left eye retina'; upside-down flip ONLY
            mappedMeridianName('temporal meridian') = 'temporal meridian';
            mappedMeridianName('nasal meridian') = 'nasal meridian';
            mappedMeridianName('superior meridian') = 'inferior meridian';
            mappedMeridianName('inferior meridian') = 'superior meridian';
        case 'right'
            %'right eye retina'; left/right &  upside-down flip
            mappedMeridianName('temporal meridian') = 'nasal meridian';
            mappedMeridianName('nasal meridian') = 'temporal meridian';
            mappedMeridianName('superior meridian') = 'inferior meridian';
            mappedMeridianName('inferior meridian') = 'superior meridian';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end
    
    eccMags = sqrt(sum(rfPositions.^2,2));
    eccAngles = atan2d(rfPositions(:,2),rfPositions(:,1));
    
    for meridianIndex = 1:numel(obj.enumeratedMeridianNames)
        [~, meridianRFDensities(meridianIndex,:)] = ...
            obj.mRGCRFSpacingAndDensityAlongMeridian(eccMags, mappedMeridianName(obj.enumeratedMeridianNames{meridianIndex}), ...
            posUnits, densityUnits,...
            'subtype', p.Results.subtype, ...
            'adjustForISETBioConeDensity',p.Results.adjustForISETBioConeDensity); 
    end
    
    % Interpolate radially
    mRGCRFDensities = obj.interpolatedValuesFromMeridianValues(meridianRFDensities, eccAngles);
    
    % Compute corresponding spacings
    mRGCSpacings = obj.spacingFromDensity(mRGCRFDensities);

end

