function importExternalConeData(obj, coneData)

    % Validate coneData struct
    validateInput(coneData);

    % Import cone positions
    switch (coneData.positionUnits)
        case 'microns'
            obj.coneRFpositionsMicrons = coneData.positions;
            if (~isempty(obj.micronsPerDegreeApproximation))
                obj.coneRFpositionsDegs = obj.coneRFpositionsMicrons/obj.micronsPerDegreeApproximation;
            else
                obj.coneRFpositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(obj.coneRFpositionsMicrons * 1e-3);
            end
        case 'degrees'
            obj.coneRFpositionsDegs = coneData.positions;
            if (~isempty(obj.micronsPerDegreeApproximation))
                obj.coneRFpositionsMicrons = obj.coneRFpositionsDegs*obj.micronsPerDegreeApproximation;
            else
                obj.coneRFpositionsMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(obj.coneRFpositionsDegs);
            end
    end
    conesNum = size(obj.coneRFpositionsMicrons,1);
    
    
    if (~isfield(coneData, 'lightGatheringApertureDiameters'))
        % Compute spacings from positions, as we do normally in cMosaic
        obj.coneRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);
        obj.coneRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsDegs);
    else
        % User supplies the cone apertures (light-collecting areas, not spacing)
        % Compute cone spacings from lightGatheringApertureDiameters
        coneApertures = reshape(coneData.lightGatheringApertureDiameters, [1 conesNum]);
        coneSpacings = coneApertures / obj.coneApertureToDiameterRatio;
        switch (coneData.positionUnits)
            case 'microns' 
                obj.coneRFspacingsMicrons = coneSpacings;
                if (~isempty(obj.micronsPerDegreeApproximation))
                    obj.coneRFspacingsDegs = obj.coneRFspacingsMicrons / obj.micronsPerDegreeApproximation;
                else
                    eccMicrons = mean(sqrt(sum(obj.coneRFpositionsMicrons.^2,2)));
                    obj.coneRFspacingsDegs = RGCmodels.Watson.convert.sizeRetinalMicronsToSizeVisualDegs(obj.coneRFspacingsMicrons, eccMicrons);
                end

            case 'degrees'
                obj.coneRFspacingsDegs = coneSpacings;
                if (~isempty(obj.micronsPerDegreeApproximation))
                    obj.coneRFspacingsMicrons = obj.coneRFspacingsDegs * obj.micronsPerDegreeApproximation; 
                else
                    eccDegs = mean(sqrt(sum(obj.coneRFpositionsDegs.^2,2)));
                    obj.coneRFspacingsMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(obj.coneRFspacingsDegs, eccDegs);
                end
        end
    end
    
    if (isfield(coneData, 'outerSegmentLengthAttenationFactors'))
        obj.importedOSLengthAttenuationFactors = reshape(coneData.outerSegmentLengthAttenationFactors, [1 conesNum]);
    end
    
    if (isfield(coneData, 'blurApertureDiameterMicrons'))
        obj.importedBlurDiameterMicrons = coneData.blurApertureDiameterMicrons;
    end
    
    
    % Compute min and max cone positions in degrees
    minEccDegs = min(obj.coneRFpositionsDegs, [], 1);
    maxEccDegs = max(obj.coneRFpositionsDegs, [], 1);
    
    % Set new eccentricity and size properties
    obj.eccentricityDegs = mean(obj.coneRFpositionsDegs,1);
    obj.sizeDegs = maxEccDegs - minEccDegs;

    % Now import cone types
    obj.coneTypes = reshape(coneData.types, [conesNum 1]);
    
    % Compute indices
    obj.lConeIndices = find(obj.coneTypes == cMosaic.LCONE_ID);
    obj.mConeIndices = find(obj.coneTypes == cMosaic.MCONE_ID);
    obj.sConeIndices = find(obj.coneTypes == cMosaic.SCONE_ID);
    obj.kConeIndices = find(obj.coneTypes == cMosaic.KCONE_ID);
    
    % Make sure all cones have been assigned an ID
    assert(conesNum==numel(obj.lConeIndices)+numel(obj.mConeIndices)+numel(obj.sConeIndices)+numel(obj.kConeIndices), ...
        'loadExternalConeData():: indices do not sum up to total cones');
    
    % Compute achieved cone densities
    achievedConeDensities = [...
        numel(obj.lConeIndices)/conesNum ...
        numel(obj.mConeIndices)/conesNum ...
        numel(obj.sConeIndices)/conesNum ...
        numel(obj.kConeIndices)/conesNum];
    
    % Update coneDensities
    obj.coneDensities = achievedConeDensities;
    
    fprintf('Achieved cone densities: L (%2.3f), M (%2.3f), S (%2.3f), K (%2.3f)\n', ...
        achievedConeDensities(1), ...
        achievedConeDensities(2), ...
        achievedConeDensities(3), ...
        achievedConeDensities(4));
end

function  validateInput(coneData)
    % Check that we have a 'positionsUnits' field
    assert(isfield(coneData, 'positionUnits'), ...
        'Expected a ''positionUnits'' field in coneData');
    assert(ismember(coneData.positionUnits, {'microns', 'degrees'}), ...
        sprintf('loadExternalConeData():: ''positionUnits'' field should be set to either ''microns'' or ''degrees'', not ''%s''.', coneData.positionUnits));    

    % Check that we have a 'positions' field
    assert(isfield(coneData, 'positions'), ...
        'loadExternalConeData():: expected a ''positions'' field in coneData');
    
    if (isfield(coneData, 'lightGatheringApertureDiameters'))
        assert(size(coneData.positions,1) == numel(coneData.lightGatheringApertureDiameters), ...
            sprintf('loadExternalConeData()::  dimensions of ''positions'' (%d) and ''lightGatheringApertureDiameters'' (%d) do not match.\n', size(coneData.positions,1), numel(coneData.lightGatheringApertureDiameters)));
    end
    
    % Check dimensionality of positions
    assert(size(coneData.positions,2) == 2, ...
         sprintf('loadExternalConeData()::  ''positions'' should be an N x 2 matrix'));
          
    % Check that we have a 'types' field
    assert(isfield(coneData, 'types'), ...
        'loadExternalConeData():: expected a ''type'' field in coneData');
     
    % Check that dimensions match
    assert(size(coneData.positions,1) == numel(coneData.types), ...
        sprintf('loadExternalConeData():: dimensions of ''positions'' (%d) and ''types'' (%d) do not match\n', size(coneData.positions,1), numel(coneData.types)));
    
end

