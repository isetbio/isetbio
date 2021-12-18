function [rfsToKeep, rfsToBeEliminated, overlappingOtherRFs] = identifyOverlappingRFs(xPos, yPos, ...
             RFpositionsMicrons, RFspacingsMicrons, maxSeparationForDeclaringOverlap)
    
    rfsNum  = size(RFpositionsMicrons,1);

    deltaTheta = 20;
    thetas = 0:deltaTheta:360-deltaTheta;
    cosThetas = cosd(thetas);
    sinThetas = sind(thetas);

    parfor iRF = 1:rfsNum-1
        r = maxSeparationForDeclaringOverlap * RFspacingsMicrons(iRF);
        xx = RFpositionsMicrons(iRF,1) + r * cosThetas;
        yy = RFpositionsMicrons(iRF,2) + r * sinThetas;
    
        overlappingRFs = [];
        for otherRF = iRF+1:rfsNum
            xCenter = RFpositionsMicrons(otherRF,1);
            yCenter = RFpositionsMicrons(otherRF,2);

            if (inpolygon(xCenter, yCenter, xx, yy))
                overlappingRFs(numel(overlappingRFs)+1) = otherRF;
            end
        end
        overlappingOtherRFs{iRF} = overlappingRFs;
    end % iRF

    problematicNodesNum = 0;
    rfsToBeEliminated = [];
    for iRF = 1:rfsNum-1
        if (~isempty(overlappingOtherRFs{iRF}))
            problematicNodesNum = problematicNodesNum + 1;
            rfsToBeEliminated = cat(2, rfsToBeEliminated , iRF);
        end
    end

    rfsToKeep = setdiff(1:rfsNum, rfsToBeEliminated);
    RFpositionsMicrons = RFpositionsMicrons(rfsToKeep,:);
    RFspacingsMicrons = RFspacingsMicrons(rfsToKeep);

    if (problematicNodesNum > 0)
        fprintf(2,'Found %d problematic nodes out of a total of %d nodes', problematicNodesNum, rfsNum);
    end

    visualizeOverlaps = false;
    if (problematicNodesNum > 0) && (visualizeOverlaps)
        hFig = figure();
        set(hFig, 'Name', sprintf('Position: %2.1f %2.f', xPos, yPos));
        maxSuplotsNum = 3;
        subplotNum = 0;
        for iRF = 1:rfsNum-1
            overlappingRFs = overlappingOtherRFs{iRF};
            if (isempty(overlappingRFs))
                continue;
            end
            
            for k = 1:numel(overlappingRFs)
                otherRF = overlappingRFs(k);
                subplotNum = subplotNum + 1;
                if (subplotNum > maxSuplotsNum^2)
                    fprintf('Skipped visualizing overlapping rf %d of %d', k, numel(overlappingRFs));
                    continue;
                end
                ax = subplot(maxSuplotsNum,maxSuplotsNum,subplotNum);
                visualizeOverlappingRFs(ax,RFspacingsMicrons, RFpositionsMicrons, cosThetas, sinThetas, iRF, otherRF);
            end % k
        end % iRF
        drawnow;
    end % if
end

function visualizeOverlappingRFs(ax, RFspacingsMicrons, RFpositionsMicrons, cosThetas, sinThetas, iRF, otherRF)
    
    r = 0.5*RFspacingsMicrons(iRF);
    xx = RFpositionsMicrons(iRF,1) + r * cosThetas;
    yy = RFpositionsMicrons(iRF,2) + r * sinThetas;
    plot(ax,xx, yy, 'k-', 'LineWidth', 1.0); hold on;
    plot(ax,RFpositionsMicrons(iRF,1), RFpositionsMicrons(iRF,2), 'ko');

    r = 0.5*RFspacingsMicrons(otherRF);
    xx = RFpositionsMicrons(otherRF,1) + r * cosThetas;
    yy = RFpositionsMicrons(otherRF,2) + r * sinThetas;
    plot(ax,xx, yy, 'r-', 'LineWidth', 1.0);
    plot(ax,RFpositionsMicrons(otherRF,1), RFpositionsMicrons(otherRF,2), 'ro');
    diffPos = RFpositionsMicrons(iRF,:) - RFpositionsMicrons(otherRF,:);    
    centerSeparation = sqrt(sum(diffPos.^2,2));
    meanSpacing = mean([RFspacingsMicrons(iRF) RFspacingsMicrons(otherRF)]);
    title(ax,sprintf('%d+%d\nseparation: %2.1fum\nmeanSpacing: %2.1fum', iRF, otherRF, centerSeparation, meanSpacing));
    axis(ax,'equal');
end
