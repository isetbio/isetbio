function reassignConeIdentities(obj, varargin)
    % Reassign the cone identities of the cone mosaic hex object
    %
    % Syntax:
    %   reassignConeIdentities(obj, [varargin])
    %
    % Description:
    %    Reassign the cone identities of the provided cone mosaic hex.
    %
    % Inputs:
    %    obj      - The cone mosaic hex object
    %
    % Outputs:
    %    None.
    %
    % Optional key/value pairs:
    %    varargin - Will contain the provided keys and values, such as...
    %    ***Needs to be filled out!***
    %
    % History:
    %    5/26/20   NPC  Fixed iRow issue, which was causing the mosaic plotting Y-coord flip 
%                   (rows grow top -> bottom, whereas Y-coords grow bottom -> top)

    p = inputParser;
    p.addParameter('sConeMinDistanceFactor', 2.5, @isnumeric);
    p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
    p.addParameter('zeroSconeDensity', false, @islogical);
    p.parse(varargin{:});

    %% Set up cone coordinates and outline
    sampledHexMosaicXaxis = obj.patternSupport(1, :, 1) + obj.center(1);
    sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2) + obj.center(2);

    
    if (p.Results.zeroSconeDensity)
        sConeIndices = find(obj.pattern == 4);
        sConeIndicesToBeReassinged = 1:numel(sConeIndices);
        [addedLconesNum, addedMconesNum] = reassignSconeIdentities(obj, ...
            sConeIndices, sConeIndicesToBeReassinged);
        sConeIndices = find(obj.pattern == 4);
        fprintf('Scones after adjustment %d: (added Lcones:%d, added Mcones:%d)\n', ...
            numel(sConeIndices), addedLconesNum, addedMconesNum);        
    else
        fprintf('Will not change the S-cone density.\n');
    end

    if (~isempty(p.Results.sConeFreeRadiusMicrons))
        sConeIndices = find(obj.pattern == 4);
        fprintf(['Scones before making an S-cone free central patch: ' ...
            '%d\n'], numel(sConeIndices));

        [iRows, iCols] = ind2sub(size(obj.pattern), sConeIndices);
        x2Microns = (1e6 * sampledHexMosaicXaxis(iCols)) .^ 2;
        y2Microns = (1e6 * sampledHexMosaicYaxis(end-iRows+1)) .^ 2;
        sConeEccentricitiesMicrons = sqrt(x2Microns(:) + y2Microns(:));
        sConeIndicesToBeReassinged = find(sConeEccentricitiesMicrons < ...
            p.Results.sConeFreeRadiusMicrons);
        [addedLconesNum, addedMconesNum] = reassignSconeIdentities(obj, ...
            sConeIndices, sConeIndicesToBeReassinged);

        sConeIndices = find(obj.pattern == 4);
        fprintf(['Scones after making an S-cone free central patch: %d' ...
            ' (added Lcones:%d, added Mcones:%d)\n'], ...
            numel(sConeIndices), addedLconesNum, addedMconesNum);
    else
        fprintf('Will not create an an S-cone free central patch.\n');
    end

    if (~isempty(p.Results.sConeMinDistanceFactor))
        keepLooping = true;
        passIndex = 1;
        while (keepLooping)
            coneLocsInMeters = [];

            sConeIndices = find(obj.pattern == 4);
            [iRows, iCols] = ind2sub(size(obj.pattern), sConeIndices);
            xxMeters = sampledHexMosaicXaxis(iCols);
            yyMeters = sampledHexMosaicYaxis(end-iRows+1);
            coneLocsInMeters(:, 1) = xxMeters;
            coneLocsInMeters(:, 2) = yyMeters;
            xxMicrons = 1e6 * xxMeters;
            yyMicrons = 1e6 * yyMeters;
            sConeCoordsMicrons = [xxMicrons(:) yyMicrons(:)];

            eccInMeters = sqrt(sum(coneLocsInMeters .^ 2, 2));
            ang = atan2(squeeze(coneLocsInMeters(:, 2)), ...
                squeeze(coneLocsInMeters(:, 1))) / pi * 180;

            if (obj.eccBasedConeDensity)
                [coneSpacing, ~, ~] = coneSizeReadData('eccentricity', ...
                    eccInMeters(:), 'angle', ang(:));
            else
                [coneSpacing, ~, ~] = coneSizeReadData('eccentricity', ...
                    0 * eccInMeters(:), 'angle', ang(:));
                if (~isempty(obj.customLambda))
                    coneSpacing = obj.customLambda / 1e6 + 0 * coneSpacing;
                end
            end

            % use 10% larger because of granularity in hex locations
            localConeSpacingMicrons = 1.1 * coneSpacing * 1e6;

            [distances, idx] = pdist2(sConeCoordsMicrons, ...
                sConeCoordsMicrons, 'euclidean', 'Smallest', 2);
            
            if (~isempty(distances)) && (size(distances, 1) == 2)
                distancesToNearestScone = squeeze(distances(2, :));
                idx = squeeze(idx(2, :));

                a = (1:size(sConeCoordsMicrons, 1));
                eliminated_a = zeros(1, size(sConeCoordsMicrons, 1));
                b = idx;

                for k = 1:numel(a)
                    % eliminate b(k) from the a-list
                    if (eliminated_a(k) == 0), eliminated_a(b(k)) = 1; end
                end

                idx = idx(eliminated_a == 0);
                distancesToNearestScone = ...
                    distancesToNearestScone(eliminated_a == 0);
                localConeSpacingMicrons = ...
                    localConeSpacingMicrons(eliminated_a == 0);
                sConeIndicesToBeReassinged = ...
                    find(distancesToNearestScone < ...
                    localConeSpacingMicrons * ...
                    p.Results.sConeMinDistanceFactor);

                [addedLconesNum, addedMconesNum] = ...
                    reassignSconeIdentities(obj, sConeIndices(idx), ...
                    sConeIndicesToBeReassinged);
                sConeIndices = find(obj.pattern == 4);
                fprintf(['Scones after pass %d: %d (added Lcones:%d, ' ...
                    'added Mcones:%d)\n'], passIndex, ...
                    numel(sConeIndices), addedLconesNum, addedMconesNum);

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

function [newLcones, newMcones] = reassignSconeIdentities(obj, ...
        sConeIndices, sConeIndicesToBeReassinged)
    % Reassign the S Cone Identities
    %
    % Syntax:
    %   reassignSconeIdentities(obj, sConeIndices, ...
    %       sConeIndicesToBeReassigned)
    %
    % Description:
    %    Reassign specific S Cone identities
    %
    % Inputs:
    %    obj                        - The cone mosaic hex object
    %    sConeIndices               - Indices of All of the S Cones
    %    sConeIndicesToBeReassigned - Indices of the relevant S Cones.
    %
    % Outputs:
    %    newLcones                  - The new L Cones
    %    newMcones                  - The new M Cones
    %
    % Optional key/value pairs:
    %    None.
    %

    LconesNum = numel(find(obj.pattern == 2));
    MconesNum = numel(find(obj.pattern == 3));
    LconesFraction = LconesNum/(LconesNum+MconesNum);
    newLcones = 0;
    newMcones = 0;
    
    randomizedIndices = randperm(numel(sConeIndicesToBeReassinged));
    
    for sIndex = 1:numel(sConeIndicesToBeReassinged) 
        
        if (LconesFraction < obj.spatialDensity(2)/(obj.spatialDensity(2)+obj.spatialDensity(3)))    
            reassignedConeIdentity = 2;
            LconesNum = LconesNum+1;
            newLcones = newLcones + 1;
        else
            reassignedConeIdentity = 3;
            MconesNum = MconesNum+1;
            newMcones = newMcones + 1;
        end
        if (reassignedConeIdentity == 2)
            if (obj.spatialDensity(2) == 0) 
                reassignedConeIdentity = 3;
                LconesNum = LconesNum-1;
                newLcones = newLcones- 1;
                MconesNum = MconesNum+1;
                newMcones = newMcones+1;
            end
        elseif (reassignedConeIdentity == 3)
            if (obj.spatialDensity(3) == 0) 
                reassignedConeIdentity = 2;
                MconesNum = MconesNum-1;
                newMcones = newMcones- 1;
                LconesNum = LconesNum+1;
                newLcones = newLcones+1;
            end
        end
        
        LconesFraction = LconesNum/(LconesNum+MconesNum);
        
        obj.pattern(sConeIndices(sConeIndicesToBeReassinged(randomizedIndices(sIndex)))) = ...
            reassignedConeIdentity;
    end
end

