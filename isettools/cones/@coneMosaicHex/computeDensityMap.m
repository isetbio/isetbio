% Method to compute the cone density of @coneMosaicHex
function [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, computeConeDensityMap)

    % Sampling interval in microns
    deltaX = 2*obj.lambdaMin*1e-6;
    
    % Sampling grid in microns
    margin = 2*deltaX;
    mosaicRadius = max([obj.width obj.height]);
    mosaicRangeX = obj.center(1) + mosaicRadius*[-1 1]  + [-obj.lambdaMin obj.lambdaMin]*1e-6;
    mosaicRangeY = obj.center(2) + mosaicRadius*[-1 1] + [-obj.lambdaMin obj.lambdaMin]*1e-6;
    gridXPos = (mosaicRangeX(1)+margin):deltaX:(mosaicRangeX(2)-margin);
    gridYPos = (mosaicRangeY(1)+margin):deltaX:(mosaicRangeY(2)-margin);
    [densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);
    
    % Allocate memory
    densityMap = zeros(numel(gridYPos), numel(gridXPos));
    
    if (strcmp(computeConeDensityMap, 'from mosaic'))
        % Determine density by computing local spacing between a cone and its 6 closest neighbors
        neigboringConesNum = 6;
        % radius of cones over which to average the local cone spacing
        averagingRadiusConesNum = 1;
        for iYpos = 1:numel(gridYPos)
            for iXpos = 1:numel(gridXPos)
                % initialize spacingInMeters
                spacingInMeters = 0;
                measuresNum = 0;
                % compute spacingInMeters over a regions of (2*averagingRadiusConesNum+1) x (2*averagingRadiusConesNum+1)
                for i = -averagingRadiusConesNum:averagingRadiusConesNum
                    ii = iXpos + i;
                    if (ii < 1) || (ii > numel(gridXPos)); continue; end
                    for j = -averagingRadiusConesNum:averagingRadiusConesNum
                        jj = iYpos + j;
                        if (jj < 1) || (jj > numel(gridYPos)); continue; end
                        % spacing between target cone and its closese neighbors
                        xoMeters = gridXPos(ii); yoMeters = gridYPos(jj);
                        nearestConeDistancesInMeters = pdist2(obj.coneLocsHexGrid, [xoMeters yoMeters], 'Euclidean', 'Smallest',neigboringConesNum);
                        % sum the mean spacing between target cone and its closese neighbors
                        spacingInMeters = spacingInMeters + mean(nearestConeDistancesInMeters);
                        measuresNum = measuresNum + 1;
                    end
                end
                % Only include measures that contain all positions within
                % the (2*averagingRadiusConesNum+1) x (2*averagingRadiusConesNum+1) area
                if (measuresNum == (2*averagingRadiusConesNum+1)^2)
                    spacingInMeters = spacingInMeters / measuresNum;
                    % spacing in mm
                    spacingInMM = spacingInMeters * 1e3;
                    % density in cones/mm2
                    densityMap(iYpos,iXpos) = (1/spacingInMM)^2;
                else
                    densityMap(iYpos,iXpos) = nan;
                end
            end
        end
    else
        [X,Y] = meshgrid(gridXPos, gridYPos);
        eccInMeters = sqrt(X.^2 + Y.^2);
        ang = atan2(Y, X)/pi*180;
        [~, ~, densities] = coneSizeReadData('eccentricity',eccInMeters(:),'angle',ang(:));
        densityMap = reshape(densities, size(X));
    end
    
end