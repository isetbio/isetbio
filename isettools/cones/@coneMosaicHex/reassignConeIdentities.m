function reassignConeIdentities(obj, varargin)

    p = inputParser;
    p.addParameter('sConeMinDistanceFactor', 3.0, @isnumeric);
    p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
    p.parse(varargin{:});
    
    
    %% Set up cone coordinates and outline
    sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
    sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);
    
    if (~isempty(p.Results.sConeFreeRadiusMicrons))
        sConeIndices = find(obj.pattern == 4);
        fprintf('Scones before making an S-cone free central patch: %d\n', numel(sConeIndices));
        
        [iRows,iCols] = ind2sub(size(obj.pattern), sConeIndices);
        x2Microns = (1e6*sampledHexMosaicXaxis(iCols)).^2;
        y2Microns = (1e6*sampledHexMosaicYaxis(iRows)).^2;
        sConeEccentricitiesMicrons = sqrt(x2Microns(:)+y2Microns(:));
        sConeIndicesToBeReassinged = find(sConeEccentricitiesMicrons < p.Results.sConeFreeRadiusMicrons);
        [addedLconesNum, addedMconesNum] = reassignSconeIdentities(obj, sConeIndices, sConeIndicesToBeReassinged);
        
        sConeIndices = find(obj.pattern == 4);
        fprintf('Scones after making an S-cone free central patch: %d (added Lcones:%d, added Mcones:%d)\n', numel(sConeIndices), addedLconesNum, addedMconesNum);
    else
        fprintf('Will not create an an S-cone free central patch.\n');
    end
    
    if (~isempty(p.Results.sConeMinDistanceFactor))
        keepLooping = true;
        passIndex = 1;
        while (keepLooping)
            coneLocsInMeters = [];
    
            sConeIndices = find(obj.pattern == 4);
            [iRows,iCols] = ind2sub(size(obj.pattern), sConeIndices);
            xxMeters = sampledHexMosaicXaxis(iCols);
            yyMeters = sampledHexMosaicYaxis(iRows);
            coneLocsInMeters(:,1) = xxMeters;
            coneLocsInMeters(:,2) = yyMeters;
            xxMicrons = 1e6*xxMeters;
            yyMicrons = 1e6*yyMeters;
            sConeCoordsMicrons = [xxMicrons(:) yyMicrons(:)];
        
            eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
            ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
            
            if (obj.eccBasedConeDensity)
                [coneSpacing, ~, ~] = coneSizeReadData(eccInMeters(:),ang(:));
            else
                [coneSpacing, ~, ~] = coneSizeReadData(0*eccInMeters(:),ang(:));
                if (~isempty(obj.customLambda))
                    coneSpacing = obj.customLambda/1e6 + 0*coneSpacing;
                end
            end
            
            % use 10% larger because of granularity in hex locations
            localConeSpacingMicrons = 1.1 * coneSpacing * 1e6;

            [distances, idx] = pdist2(sConeCoordsMicrons, sConeCoordsMicrons, 'euclidean', 'Smallest', 2);
            
            if (~isempty(distances)) && (size(distances,1) == 2)
                distancesToNearestScone = squeeze(distances(2,:));
                idx = squeeze(idx(2,:));

                a = (1:size(sConeCoordsMicrons,1));
                eliminated_a = zeros(1,size(sConeCoordsMicrons,1));
                b = idx;

                for k = 1:numel(a)
                    % eliminate b(k) from the a-list
                    if (eliminated_a(k) == 0)
                        eliminated_a(b(k)) = 1;
                    end
                end
        
                idx = idx(eliminated_a == 0);
                distancesToNearestScone = distancesToNearestScone(eliminated_a == 0);
                localConeSpacingMicrons = localConeSpacingMicrons(eliminated_a == 0);
                sConeIndicesToBeReassinged = find(distancesToNearestScone < localConeSpacingMicrons*p.Results.sConeMinDistanceFactor);

                [addedLconesNum, addedMconesNum] = reassignSconeIdentities(obj, sConeIndices(idx), sConeIndicesToBeReassinged);
                sConeIndices = find(obj.pattern == 4);
                fprintf('Scones after pass %d: %d (added Lcones:%d, added Mcones:%d)\n', passIndex, numel(sConeIndices), addedLconesNum, addedMconesNum);
                
                if ((addedLconesNum == 0) && (addedMconesNum == 0))
                    keepLooping = false;
                else
                	passIndex = passIndex + 1;
                end
            else
                keepLooping = false;
            end 
        end % keepLooping
    else
        fprintf('Will not adjust the minimum S-cone separation\n.');
    end % p.Results.sConeMinDistanceFactor
    
end

function [newLcones, newMcones] = reassignSconeIdentities(obj, sConeIndices, sConeIndicesToBeReassinged)
    newLcones = 0;
    newMcones = 0;
    theRandomNumbers = rand(1,numel(sConeIndicesToBeReassinged));
    MconeFraction = obj.spatialDensity(3)/(obj.spatialDensity(2)+obj.spatialDensity(3));
        
    for sIndex = 1:numel(sConeIndicesToBeReassinged) 
        if (theRandomNumbers(sIndex) < MconeFraction)
            reassignedConeIdentity = 3;
            newMcones = newMcones + 1;
        else
            reassignedConeIdentity = 2;
            newLcones = newLcones + 1;
        end
        obj.pattern(sConeIndices(sConeIndicesToBeReassinged)) = reassignedConeIdentity;
    end
end
