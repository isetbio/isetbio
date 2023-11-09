function customConeFundamentals = coneFundamentalsAtTargetPositionWithinConeMosaic(...
    theConeMosaic, theOptics, targetRegionPositionDegs, targetRegionSizeDegs, maxConesNumForAveraging)

    % Find cone indices within the target region
    if (numel(targetRegionSizeDegs) == 2)
        widthDegs = targetRegionSizeDegs(1);
        heightDegs = targetRegionSizeDegs(2);
    else
        widthDegs = targetRegionSizeDegs(1);
        heightDegs = targetRegionSizeDegs(1);
    end

    theStimulusRegion = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', targetRegionPositionDegs, ...
            'width', widthDegs, ...
            'height', heightDegs, ...
            'rotation', 0.0...
        ));
    targetConeIndices = theStimulusRegion.indicesOfPointsInside(theConeMosaic.coneRFpositionsDegs);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.LCONE_ID);
    indicesOfLconesWithinTargetRegion = targetConeIndices(idx);
    dd = theConeMosaic.coneRFpositionsDegs(indicesOfLconesWithinTargetRegion,:);
    dd = sqrt(sum((bsxfun(@minus, dd, targetRegionPositionDegs)).^2,2));
    [~,iidx] = sort(dd, 'ascend');
    indicesOfLconesWithinTargetRegion = targetConeIndices(idx(iidx));
    indicesOfLconesWithinTargetRegion = indicesOfLconesWithinTargetRegion(1:maxConesNumForAveraging);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.MCONE_ID);
    indicesOfMconesWithinTargetRegion  = targetConeIndices(idx);
    dd = theConeMosaic.coneRFpositionsDegs(indicesOfMconesWithinTargetRegion,:);
    dd = sqrt(sum((bsxfun(@minus, dd, targetRegionPositionDegs)).^2,2));
    [~,iidx]  = sort(dd, 'ascend');
    indicesOfMconesWithinTargetRegion = targetConeIndices(idx(iidx));
    indicesOfMconesWithinTargetRegion = indicesOfMconesWithinTargetRegion(1:maxConesNumForAveraging);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.SCONE_ID);
    if (~isempty(idx))
        indicesOfSconesWithinTargetRegion  = targetConeIndices(idx);
        dd = theConeMosaic.coneRFpositionsDegs(indicesOfSconesWithinTargetRegion,:);
        dd = sqrt(sum((bsxfun(@minus, dd, targetRegionPositionDegs)).^2,2));
        [~,iidx]  = sort(dd, 'ascend');
        indicesOfSconesWithinTargetRegion = targetConeIndices(idx(iidx));
        indicesOfSconesWithinTargetRegion = indicesOfSconesWithinTargetRegion(1:maxConesNumForAveraging);
    else
        indicesOfSconesWithinTargetRegion = [];
        fprintf(2, 'No S-cones in the region around (%2.2f,%2.2f) degs. Will use the SS-2 S-cone fundamental.\n', targetRegionPositionDegs(1), targetRegionPositionDegs(2))
    end

    fprintf('Target region contains %d L-cones, %d M-cones and %d S-cones\n', ...
        numel(indicesOfLconesWithinTargetRegion), numel(indicesOfMconesWithinTargetRegion), numel(indicesOfSconesWithinTargetRegion));

    % Initialize
    customConeFundamentals.wavelengthSupport = oiGet(theOptics, 'wave');
    customConeFundamentals.quantalExcitationSpectra = zeros(numel(customConeFundamentals.wavelengthSupport),3);
    
    StockmanSharpe2DegConeFundamentals = ...
        ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), customConeFundamentals.wavelengthSupport);

    % Compute the cone mosaic's response to a series of monochromatic, 
    % spatially-uniform images with constant power in photons/sec-nm
    pixelsNum = 256;
    
    for iMonoChromaticBand = 1:numel(customConeFundamentals.wavelengthSupport)
        % Set up a spatially uniform dummy scene
        % (black body spectrum of 5000 degK).
        scene  = sceneCreate('uniform bb', pixelsNum, 5000, customConeFundamentals.wavelengthSupport);
    
        % Set the scene FOV to 1 deg
        fovDegs = 1;
        scene = sceneSet(scene,'fov', fovDegs);

        % Rewrite the scene with the desired monochromatic spectrum
        % allowing photons only at the current monomchromatic wavelength
        energyWithinBand = 1;
        energy = sceneGet(scene,'energy')*0;
        energy(:,:,iMonoChromaticBand) = energyWithinBand * ones(pixelsNum,pixelsNum);
        scene = sceneSet(scene,'energy',energy);

        % Compute the retinal image of the monochromatic scene
        theMonochromaticOpticalImage = oiCompute(scene, theOptics);

        % Compute the L-, M-, and S-cone exciations to the monochromatic scene
        theConeExcitations = theConeMosaic.compute(...
            theMonochromaticOpticalImage, ...
            'opticalImagePositionDegs', targetRegionPositionDegs, ...
            'nTrials', 1);

        theConeExcitations = theConeExcitations(:);
        
        % Compute the quantal L-cone excitations at this wavelength
        customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.LCONE_ID) = ...
            mean(theConeExcitations(indicesOfLconesWithinTargetRegion)) / energyWithinBand;

        % Compute the quantal M-cone excitatiomns at this wavelength
        customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.MCONE_ID) = ...
            mean(theConeExcitations(indicesOfMconesWithinTargetRegion)) / energyWithinBand;

        % Compute the quantal S-cone excitations at this wavelength
        if (isempty(indicesOfSconesWithinTargetRegion))
            % Just use the SS-2 S-cone fundamental since there are no
            % S-cones on the analyzed path
            customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.SCONE_ID) = ...
                StockmanSharpe2DegConeFundamentals(iMonoChromaticBand, cMosaic.SCONE_ID);
        else
            customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.SCONE_ID) = ...
                mean(theConeExcitations(indicesOfSconesWithinTargetRegion)) / energyWithinBand;
        end

    end % iMonochromaticBand 

    figure(99); clf;
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,1), 'r-');
    hold on;
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,2), 'g-');
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,3), 'b-');
    drawnow;

    for iCone = 1:3
        customConeFundamentals.quantalExcitationSpectra(:,iCone) = ...
        customConeFundamentals.quantalExcitationSpectra(:,iCone) / max(squeeze(customConeFundamentals.quantalExcitationSpectra(:,iCone)));
    end

end