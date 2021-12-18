function testMosaics

    % The higher this overlap the more false alarms
    maxSeparationForDeclaringOverlap = 0.4;

    doEntireMosaic = ~true;
    if (doEntireMosaic)
        whichEye = 'right eye';
        sourceLatticeSizeDegs = 58;
        p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', whichEye);
        theMosaicFileName = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
        fprintf('Importing data from %s\n', theMosaicFileName);
        load(theMosaicFileName, 'fovDegs', 'neuronType', 'params', 'whichEye', 'rfPositions');
        fprintf('Computing spacings for %d cones', size(rfPositions,1));
        rfSpacings = RGCmodels.Watson.convert.positionsToSpacings(rfPositions);
        [rfsToKeep, rfsToBeEliminated, overlappingOtherRFs] = cMosaic.identifyOverlappingRFs(0,0, rfPositions, rfSpacings, maxSeparationForDeclaringOverlap);

        % Replace the position/spacing of the other overlapping RF with the
        % average position/spacing of the iRF and the otherRF
        rfsNum = size(rfPositions,1);
        for iRF = 1:rfsNum-1
             otherRFs = overlappingOtherRFs{iRF};
             if (~isempty(otherRFs))
                otherRF = otherRFs(1);
                rfPositions(otherRF,:) = 0.5*(rfPositions(otherRF,:) + rfPositions(iRF,:));
             end
        end
        
        fprintf('Saving %d out of %d cone positions to %s\n', numel(rfsToKeep), size(rfPositions,1),theMosaicFileName);
        rfPositions = rfPositions(rfsToKeep,:);
        save(theMosaicFileName, 'fovDegs', 'neuronType', 'params', 'whichEye', 'rfPositions');
    else
        sizeDegs = [2 2];
        for xPos = -10:1:10
            for yPos = -10:1:10
                c = cMosaic('eccentricityDegs', [xPos yPos], 'sizeDegs', sizeDegs);
                cMosaic.identifyOverlappingRFs(xPos, yPos, c.coneRFpositionsMicrons, c.coneRFspacingsMicrons, maxSeparationForDeclaringOverlap);
            end
        end
    end

end

