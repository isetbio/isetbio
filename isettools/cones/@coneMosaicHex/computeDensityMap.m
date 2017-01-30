% Method to compute the cone density of @coneMosaicHex
function [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, computeConeDensityMap)

    deltaX = 3*obj.lambdaMin*1e-6;
    areaHalfWidth = deltaX*4;

    mosaicRangeX = obj.center(1) + obj.width/2*[-1 1]  + [-obj.lambdaMin obj.lambdaMin]*1e-6;
    mosaicRangeY = obj.center(2) + obj.height/2*[-1 1] + [-obj.lambdaMin obj.lambdaMin]*1e-6;

    gridXPos = (mosaicRangeX(1)+areaHalfWidth):deltaX:(mosaicRangeX(2)-areaHalfWidth);
    gridYPos = (mosaicRangeY(1)+areaHalfWidth):deltaX:(mosaicRangeY(2)-areaHalfWidth);
    measurementAreaInMM2 = (2*areaHalfWidth*1e3)^2;
    densityMap = zeros(numel(gridYPos), numel(gridXPos));

    if (strcmp(computeConeDensityMap, 'from mosaic'))
        for iYpos = 1:numel(gridYPos)
            for iXpos = 1:numel(gridXPos)
                xo = gridXPos(iXpos);
                yo = gridYPos(iYpos);
                conesWithin = numel(find( ...
                    obj.coneLocsHexGrid(:,1) >= xo-areaHalfWidth  & ...
                    obj.coneLocsHexGrid(:,1) <= xo+areaHalfWidth  & ...
                    obj.coneLocsHexGrid(:,2) >= yo-areaHalfWidth  & ...
                    obj.coneLocsHexGrid(:,2) <= yo+areaHalfWidth ));
                densityMap(iYpos,iXpos) = conesWithin / measurementAreaInMM2;
            end
        end
    else
        [X,Y] = meshgrid(gridXPos, gridYPos);
        eccInMeters = sqrt(X.^2 + Y.^2);
        ang = atan2(Y, X)/pi*180;
        [~, ~, densities] = coneSize(eccInMeters(:),ang(:));
        densityMap = reshape(densities, size(X));
    end

    gridXPos = mosaicRangeX(1) + (gridXPos - min(gridXPos))/(max(gridXPos) - min(gridXPos))*(mosaicRangeX(2)-mosaicRangeX(1));
    gridYPos = mosaicRangeY(1) + (gridYPos - min(gridYPos))/(max(gridYPos) - min(gridYPos))*(mosaicRangeY(2)-mosaicRangeY(1));
    [densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);
end