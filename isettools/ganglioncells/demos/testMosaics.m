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
        sizeDegs = [4 4];
        for xPos = 0 %-16:2:16
            for yPos = 0% -16:2:16
                fprintf('Testing mosaic at %2.0f %2.1f\n', xPos, yPos);
                c = cMosaic('eccentricityDegs', [xPos yPos], 'sizeDegs', sizeDegs);
                [~, rfsToBeEliminated, ~] = cMosaic.identifyOverlappingRFs(...
                    xPos, yPos, c.coneRFpositionsMicrons, c.coneRFspacingsMicrons, ...
                    maxSeparationForDeclaringOverlap);
                fprintf('\tPost initialization test: Found %2.0f overlapping elements\n\n', numel(rfsToBeEliminated));
                if (~isempty(rfsToBeEliminated))
                    c.visualize('outlinedConesWithIndices', rfsToBeEliminated);
                    fprintf('Found %d overlapping cones.\n', numel(rfsToBeEliminated));
                end
            end
        end
    end

end

