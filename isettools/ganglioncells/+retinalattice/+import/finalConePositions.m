function rfPositionsMicrons = finalConePositions(sourceLatticeSizeDegs, eccDegs, sizeDegs, whichEye)

    % Convert degs to retinal microns
    eccMicrons  = 1000*RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));

    % Load final cone positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', whichEye);
    theMosaicFileName = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
    load(theMosaicFileName, 'rfPositions');
    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));

    eliminateOvelappingElements = true;
    if (eliminateOvelappingElements)
        for pass = 1:3
            fprintf('Checking for overlapping elements within a population of %2.0f elements (PASS #%d)...\n', ...
                size(rfPositionsMicrons,1), pass);
            tic
    
            % Compute spacings
            rfSpacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(rfPositionsMicrons);
    
            % Check for overlapping elements
            maxSeparationForDeclaringOverlap = 0.4;
            
            % Identify elements that overlap
            [rfsToKeep, rfsToBeEliminated, overlapingRFindex] = cMosaic.identifyOverlappingRFs(0,0, ...
                rfPositionsMicrons, rfSpacingsMicrons, maxSeparationForDeclaringOverlap);
    
            % Replace the position/spacing of the other overlapping RF with the
            % average position/spacing of the iRF and the otherRF
            for k = 1:numel(rfsToBeEliminated)
                iRF = rfsToBeEliminated(k);
                pos1 = rfPositionsMicrons(iRF,:);
                overlappingRF = overlapingRFindex(iRF);
                pos2 = rfPositionsMicrons(overlappingRF,:);
                rfPositionsMicrons(overlappingRF,:) = 0.5*(pos1+pos2);
            end
        
            % Only keep the rfs that are non-overlapping
            rfPositionsMicrons = rfPositionsMicrons(rfsToKeep,:);
            fprintf('Overlapping element detection  took %2.1f seconds\n',toc);
            if (numel(rfsToBeEliminated)>0)
                fprintf(2, 'Eliminated %2.0f overlapping elements during PASS #%d.\n', ...
                    numel(rfsToBeEliminated), pass);
            end
        end % pass

end

