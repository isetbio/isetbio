function testMosaics

    % The higher this overlap the more false alarms
    maxSeparationForDeclaringOverlap = 0.5;

    doEntireMosaic = ~true;
    if (doEntireMosaic)
        whichEye = 'right eye';
        sourceLatticeSizeDegs = 58;
        p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', whichEye);
        theMosaicFileName = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
        fprintf('Importing data from %s\n', theMosaicFileName);
        load(theMosaicFileName, 'fovDegs', 'neuronType', 'params', 'whichEye', 'rfPositions');
        fprintf('Computing spacings for %d cones', size(rfPositions,1));

        % Initialize
        maxPassesNum = 4; pass = 0;
        rfsToKeep = []; previousRFsNum = size(rfPositions,1);
        
        while (pass < maxPassesNum) && (numel(rfsToKeep) < previousRFsNum)
            pass = pass + 1;
            fprintf('Checking for overlapping elements within a population of %2.0f elements (PASS #%d)...\n', ...
                size(rfPositions,1), pass);
            tic
    
            % Compute spacings
            rfSpacings = RGCmodels.Watson.convert.positionsToSpacings(rfPositions);
            
            % Identify elements that overlap
            [rfsToKeep, rfsToBeEliminated, overlapingRFindex] = cMosaic.identifyOverlappingRFs(0,0, ...
                rfPositions, rfSpacings, maxSeparationForDeclaringOverlap);
    
            % Replace the position/spacing of the other overlapping RF with the
            % average position/spacing of the iRF and the otherRF
            for k = 1:numel(rfsToBeEliminated)
                iRF = rfsToBeEliminated(k);
                pos1 = rfPositions(iRF,:);
                overlappingRF = overlapingRFindex(iRF);
                pos2 = rfPositions(overlappingRF,:);
                rfPositions(overlappingRF,:) = 0.5*(pos1+pos2);
            end
        
            % Only keep the rfs that are non-overlapping
            rfPositions = rfPositions(rfsToKeep,:);
            fprintf('Overlapping element detection  took %2.1f seconds\n',toc);
            if (numel(rfsToBeEliminated)>0)
                fprintf(2, 'Eliminated %2.0f overlapping elements during PASS #%d.\n', ...
                    numel(rfsToBeEliminated), pass);
            end
        end % pass

        fprintf('Saving non-overlapping data out to %s\n', theMosaicFileName);
        save(theMosaicFileName, 'fovDegs', 'neuronType', 'params', 'whichEye', 'rfPositions');
    else
        sizeDegs = 4*[1 1];
        for xPos = -15:3:15
            for yPos = -15:3:15
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

