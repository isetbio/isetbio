function rfPositionsMicrons = finalConePositions(sourceLatticeSizeDegs, eccDegs, sizeDegs, whichEye, overlappingConeFractionForElimination)

    % Convert degs to retinal microns
    eccMicrons  = 1000*RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));

    % Load final cone positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', whichEye);
    theMosaicFileName = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
    fprintf('Loading cone mosaic data from %s.\n', theMosaicFileName);
    load(theMosaicFileName, 'rfPositions');

    % Reverse the polarity
    rfPositions = -rfPositions;

    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));

    if (~isempty(overlappingConeFractionForElimination))
        % Check for overlapping elements within this max separation
        maxSeparationForDeclaringOverlap = overlappingConeFractionForElimination;
        fprintf('Will check and eliminate overlapping cones (threshold: %2.2f)\n', maxSeparationForDeclaringOverlap);
           
        % Initialize
        maxPassesNum = 4; pass = 0;
        rfsToKeep = []; previousRFsNum = size(rfPositionsMicrons,1);
        
        while (pass < maxPassesNum) && (numel(rfsToKeep) < previousRFsNum)
            
        	pass = pass + 1;
            fprintf('Checking for overlapping elements within a population of %2.0f elements (PASS #%d)...\n', ...
                size(rfPositionsMicrons,1), pass);
            tic
    
            % Compute spacings
            rfSpacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(rfPositionsMicrons);
            
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
            previousRFsNum = size(rfPositionsMicrons,1);
            rfPositionsMicrons = rfPositionsMicrons(rfsToKeep,:);
            fprintf('Overlapping element detection  took %2.1f seconds\n',toc);
            if (numel(rfsToBeEliminated)>0)
                fprintf(2, 'Eliminated %2.0f overlapping elements during PASS #%d.\n', ...
                    numel(rfsToBeEliminated), pass);
            end

        end % pass
    end

end

