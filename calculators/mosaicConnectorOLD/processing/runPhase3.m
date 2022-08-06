% Phase 3: Assign cone types to the coneMosaic
function runPhase3(runParams)

    % Load data
    load(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.inputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios');
    
    % Assign types (L,M,S) in the cone mosaic
    coneTypes = assignConeTypes(conePositionsMicrons, coneSpacingsMicrons, ...
            runParams.tritanopicAreaDiameterMicrons, runParams.relativeSconeSpacing, runParams.LtoMratio);

     % Instantiate a plotlab object
    plotlabOBJ = plotlab();

    % Apply the default plotlab recipe overriding 
    % the color order and the figure size
    figHeightInches = 15;
    plotlabOBJ.applyRecipe(...
        'renderer', 'painters', ... %'opengl', ...
        'axesBox', 'on', ...
        'colorOrder', [0 0 0; 1 0 0.5], ...
        'axesTickLength', [0.015 0.01]/2,...
        'axesFontSize', 22, ...
        'figureWidthInches', figHeightInches*2, ...
        'figureHeightInches', figHeightInches);
    
    % Visualize central 2x1 degs
    roi = struct('xo', 0.0, 'yo', 0.0, 'width', 2, 'height', 1);
    visualizeConeMosaic(conePositionsMicrons, coneTypes, roi, plotlabOBJ)
     
    % Save mosaic data within the roiRadiusDeg
    save(fullfile(runParams.outputDir, sprintf('%s.mat',runParams.outputFile)), ...
            'conePositionsMicrons', 'coneSpacingsMicrons', 'coneTypes', ...
            'RGCRFPositionsMicrons', 'RGCRFSpacingsMicrons', ...
            'desiredConesToRGCratios');
end

function coneTypes = assignConeTypes(conePositionsMicrons, coneSpacingsMicrons, ...
            tritanopicAreaDiameterMicrons, relativeSconeSpacing, LtoMratio)
   
    conesNum = size(conePositionsMicrons,1);
       
    % Reserve all cones within the tritanopic area to be either L or M
    ecc = sqrt(sum(conePositionsMicrons.^2,2));
    fovealLorMconeIndices = find(ecc <= 0.5*tritanopicAreaDiameterMicrons);
    
    % Reserve cones outside the tritanopic are to be either L or M, except
    % cones that are separated by relativeSconeSpacing*local cone spacing,
    % which will remain S
    peripheralConeIndices = setdiff(1:conesNum, fovealLorMconeIndices);
    idx = determineLMconeIndices(conePositionsMicrons(peripheralConeIndices,:), coneSpacingsMicrons(peripheralConeIndices), relativeSconeSpacing);
    peripheralLMConeIndices = peripheralConeIndices(idx);
    
    % Assign type to all peripheralConeIndices (to L for now)
    LMconeIndices = [fovealLorMconeIndices(:); peripheralLMConeIndices(:)];
    
    % Assing L and M-cone indices
    p = rand(1,numel(LMconeIndices));
    if (isinf(LtoMratio))
        idx = 1:numel(LMconeIndices);
    else
        idx = find(p<LtoMratio/(1+LtoMratio));
    end
    LconeIndices = LMconeIndices(idx);
    MconeIndices = setdiff(LMconeIndices, LconeIndices);
    SconeIndices = setdiff(1:conesNum, LMconeIndices);
    assert(conesNum-numel(LconeIndices)-numel(MconeIndices)-numel(SconeIndices)==0, ...
        'Indices do not sum up to total cones');
    
    % Return cone types
    coneTypes = zeros(1,conesNum);
    coneTypes(LconeIndices) = 2;
    coneTypes(MconeIndices) = 3;
    coneTypes(SconeIndices) = 4;
    
end


function LMconeIndices = determineLMconeIndices(conePositionsMicrons, coneSpacingsMicrons, relativeSconeSpacing)
    conesNum = size(conePositionsMicrons,1);
    
    % Compute ecc of all cones
    ecc = sqrt(sum(conePositionsMicrons.^2,2));
    
    % Compute distances between each cone and its closest 100 cones
    [d, i] = pdist2(conePositionsMicrons, conePositionsMicrons, 'euclidean', 'smallest', 200);
    
    % Remove the distance to the cone itself
    d = d(2:end,:);
    i = i(2:end,:);
    
    % Go through all cones assigning as S-cones those that are no closer
    % than coneSpacingsMicrons(coneIndex)*relativeSconeSpacing from each other
    coneIndex = 1;
    LMconeIndices = [];
    SconeIndices = [];
    remainingConeIndices = 1:conesNum;
    
    while (numel(remainingConeIndices)>0)
        % Leave the type of the current coneIndex as S.
        SconeIndices = cat(2, SconeIndices, coneIndex);
        % This means all cones around it within the exclusion radius must be non S
        currentExclusionRadius = coneSpacingsMicrons(coneIndex)*relativeSconeSpacing;
        distancesToNearbyCones = d(:,coneIndex);
        idx = find(distancesToNearbyCones < currentExclusionRadius);
        LMconeIndices = cat(1, LMconeIndices, squeeze(i(idx, coneIndex)));
        % Keep a count of the remaining cone indices that need to be visited
        remainingConeIndices = setdiff(remainingConeIndices, SconeIndices);
        remainingConeIndices = setdiff(remainingConeIndices, LMconeIndices);
        % Next cone to visit
        [~,idx] = min(ecc(remainingConeIndices));
        coneIndex = remainingConeIndices(idx);
    end
    LMconeIndices = unique(LMconeIndices);
end
