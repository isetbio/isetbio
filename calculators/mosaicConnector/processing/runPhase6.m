% Phase 6: Compute cone connections to mRGC  RF surrounds
function runPhase6(runParams)

    % Load data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios', 'midgetRGCconnectionMatrix');
     
 
    % Compute RF center radius for each mRGC
    rfCenterRadii = 
    rfEccDegs = 
    % Call 
    ck = CronerKaplanRGCModel();
    
    modelCenterRadii = ck.centerRadiusFunction(obj.centerRadiusParams, eccDegs);
    modelSurroundRadii = obj.surroundRadiusFunction(obj.surroundRadiusParams, eccDegs);
    modelCenterSurroundRadiusRatios = centerRadii ./ surroundRadii;
    
    for rfIndex = 1:numel(rfCenterRadii)
        rfCenterRadius = rfCenterRadii(rfIndex)
        
    end
    
end