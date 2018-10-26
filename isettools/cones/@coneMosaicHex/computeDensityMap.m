function [densityMap, densityMapSupportX, densityMapSupportY] = ...
    computeDensityMap(obj, computeConeDensityMap)


% Sampling interval in microns
deltaX = 3 * obj.lambdaMin * 1e-6;
margin = 0 * obj.lambdaMin * 1e-6;

mosaicSize = 0.5*max([obj.width obj.height]);
gridXPos = 0:deltaX:(mosaicSize-margin);
gridXPos = cat(2,-fliplr(gridXPos), gridXPos(2:end));
gridYPos = gridXPos;
[densityMapSupportX,densityMapSupportY] = meshgrid(gridXPos, gridYPos);

if (strcmp(computeConeDensityMap, 'from mosaic'))
    % figure(111);
    % clf;
    % plot(obj.coneLocsHexGrid(:,1), obj.coneLocsHexGrid(:,2), 'ko');
    % hold on;
    % plot(densityMapSupportX(:), densityMapSupportY(:), 'r*');
    % axis 'square'


    meanDistanceInMM = zeros(1,numel(densityMapSupportX));

    for iPos = 1:numel(densityMapSupportX)
        % Find the cone closest to the grid point
        gridPos = [densityMapSupportX(iPos) densityMapSupportY(iPos)];
        distances = sum(bsxfun(@minus, obj.coneLocsHexGrid, gridPos).^2,2);
        [~,targetConeIndex] = min(distances);

        % All the conex minus the target cone
        tmpConeLocs = obj.coneLocsHexGrid;
        tmpConeLocs(targetConeIndex,:) = nan;
        
        % Find the mean distances of the target cone to its 6 neighboring cones
        neigboringConesNum = 4;
        nearestConeDistancesInMeters = pdist2(tmpConeLocs, ...
            obj.coneLocsHexGrid(targetConeIndex,:), ...
            'Euclidean', 'Smallest', neigboringConesNum);

        meanDistanceInMM(1,iPos) = mean(nearestConeDistancesInMeters) * 1000;
    end

    fprintf('Min distance: %2.2f microns\n', min(meanDistanceInMM)*1000);
    
    densityMap = (1./meanDistanceInMM).^2;
    densityMap = reshape(densityMap, [size(densityMapSupportX,1) size(densityMapSupportX,2)]);

    hsize = 3;
    sigma = 1.6/2;
    smoothingKernel = fspecial('gaussian',hsize,sigma);
    densityMapSmoothed = conv2(densityMap, smoothingKernel, 'same');
    densityMap = densityMapSmoothed;
    
%     figure(221);
%     imagesc(smoothingKernel)
%         axis 'image'
%     colormap(gray)
% 
%     figure(222); clf;
%     subplot(1,2,1)
%     imagesc(densityMap);
%     axis 'image';
%     colorbar
%     subplot(1,2,2)
%     imagesc(densityMapSmoothed);
%     axis 'image';
%     colormap(gray);
%     colorbar
%     pause
    
else
    eccInMeters = sqrt(densityMapSupportX .^ 2 + densityMapSupportY .^ 2);
    ang = atan2(densityMapSupportY, densityMapSupportX) / pi * 180;
    [~, ~, densities] = coneSizeReadData(...
        'eccentricity', eccInMeters(:), 'angle', ang(:));
    densityMap = reshape(densities, size(densityMapSupportX));
end
end


function [densityMap, densityMapSupportX, densityMapSupportY] = ...
    computeDensityMapOLD(obj, computeConeDensityMap)
% Method to compute the cone density of @coneMosaicHex
%
% Syntax:
%   [densityMap, densityMapSupportX, densityMapSupportY] = ...
%       computeDensityMap(obj, computeConeDensityMap)
%
% Description:
%    This is the method to compute the cone density of a cone mosaic hex
%
% Inputs:
%    obj                   - A cone mosaic hex object
%    computeConeDensityMap - The cone density map
%
% Outputs:
%    densityMap            - The created density map
%    densitySupportX       - The X-axis support
%    densitySupportY       - The Y-axis support
%
% Optional key/value pairs:
%    None.
%


% Sampling interval in microns
deltaX = 2 * obj.lambdaMin * 1e-6;

% Sampling grid in microns
margin = 2 * deltaX;
mosaicRadius = max([obj.width obj.height]);
mosaicRangeX = obj.center(1) + mosaicRadius * [-1 1] + ...
    [-obj.lambdaMin obj.lambdaMin] * 1e-6;
mosaicRangeY = obj.center(2) + mosaicRadius * [-1 1] + ...
    [-obj.lambdaMin obj.lambdaMin] * 1e-6;
gridXPos = (mosaicRangeX(1) + margin):deltaX:(mosaicRangeX(2) - margin);
gridYPos = (mosaicRangeY(1) + margin):deltaX:(mosaicRangeY(2) - margin);
[densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);

% Allocate memory
densityMap = zeros(numel(gridYPos), numel(gridXPos));

if (strcmp(computeConeDensityMap, 'from mosaic'))
    % Determine density by computing local spacing between a cone and its
    % closest neighbors
    neigboringConesNum = 6;
    % radius of cones over which to average the local cone spacing
    averagingRadiusConesNum = 1;
    for iYpos = 1:numel(gridYPos)
        for iXpos = 1:numel(gridXPos)
            % initialize spacingInMeters
            spacingInMeters = 0;
            measuresNum = 0;
            % compute spacingInMeters over the regions of:
            %   (2 * averagingRadiusConesNum + 1) x ...
            %   (2 * averagingRadiusConesNum + 1)
            for i = -averagingRadiusConesNum:averagingRadiusConesNum
                ii = iXpos + i;
                if (ii < 1) || (ii > numel(gridXPos)); continue; end
                for j = -averagingRadiusConesNum:averagingRadiusConesNum
                    jj = iYpos + j;
                    if (jj < 1) || (jj > numel(gridYPos)); continue; end
                    % spacing between target cone and its closest neighbors
                    xoMeters = gridXPos(ii); yoMeters = gridYPos(jj);
                    nearestConeDistancesInMeters = pdist2(...
                        obj.coneLocsHexGrid, [xoMeters yoMeters], ...
                        'Euclidean', 'Smallest', neigboringConesNum);
                    % sum the mean spacing between target cone and its
                    % closest neighbors
                    spacingInMeters = spacingInMeters + ...
                        mean(nearestConeDistancesInMeters);
                    measuresNum = measuresNum + 1;
                end
            end
            % Only include measures that contain all positions within the
            % (2 * averagingRadiusConesNum + 1) x ...
            %    (2 * averagingRadiusConesNum + 1) area
            if (measuresNum == (2 * averagingRadiusConesNum + 1) ^ 2)
                spacingInMeters = spacingInMeters / measuresNum;
                % spacing in mm
                spacingInMM = spacingInMeters * 1e3;
                % density in cones/mm2
                densityMap(iYpos, iXpos) = (1 / spacingInMM) ^ 2;
            else
                densityMap(iYpos, iXpos) = nan;
            end
        end
    end
else
    [X, Y] = meshgrid(gridXPos, gridYPos);
    eccInMeters = sqrt(X .^ 2 + Y .^ 2);
    ang = atan2(Y, X) / pi * 180;
    [~, ~, densities] = coneSizeReadData(...
        'eccentricity', eccInMeters(:), 'angle', ang(:));
    densityMap = reshape(densities, size(X));
end

end