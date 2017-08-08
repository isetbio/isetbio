% Method to compute the cone density of @coneMosaicHex
function [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, computeConeDensityMap)

    deltaX = 3*obj.lambdaMin*1e-6;
    margin = 4*deltaX;
    mosaicRangeX = obj.center(1) + obj.width*[-1 1]  + [-obj.lambdaMin obj.lambdaMin]*1e-6;
    mosaicRangeY = obj.center(2) + obj.height*[-1 1] + [-obj.lambdaMin obj.lambdaMin]*1e-6;
    gridXPos = (mosaicRangeX(1)+margin):deltaX:(mosaicRangeX(2)-margin);
    gridYPos = (mosaicRangeY(1)+margin):deltaX:(mosaicRangeY(2)-margin);
    densityMap = zeros(numel(gridYPos), numel(gridXPos));

    halfConeApertureMicrons = 0.5*obj.lambdaMin;
    coneLocsHexGrid = obj.coneLocsHexGrid * 1e6;
    coneLocsHexGridX = squeeze(coneLocsHexGrid(:,1));
    coneLocsHexGridY = squeeze(coneLocsHexGrid(:,2));
    
    if (strcmp(computeConeDensityMap, 'from mosaic'))
        for iYpos = 1:numel(gridYPos)
            for iXpos = 1:numel(gridXPos)
                xoMeters = gridXPos(iXpos);
                yoMeters = gridYPos(iYpos);
                xoMicrons = xoMeters * 1e6;
                yoMicrons = yoMeters * 1e6;
                eccentricityInMeters = sqrt(xoMeters^2 + yoMeters^2);
                ang = atan2(yoMeters, xoMeters)/pi*180;
                [coneSpacingInMeters aperture density] = coneSize(eccentricityInMeters,ang);
                % count cones within a region 12x12 times the local cone spacing
                % so our counting area grows proprotionally with the local cone spacing
                measurementAreaHalfWidthMicrons = 6*coneSpacingInMeters*1e6;
                measurementAreaInMM2 = (2*measurementAreaHalfWidthMicrons*1e-3)^2;
                conesWithin = numel(find( ...
                    coneLocsHexGridX-halfConeApertureMicrons >= xoMicrons-measurementAreaHalfWidthMicrons  & ...
                    coneLocsHexGridX+halfConeApertureMicrons <= xoMicrons+measurementAreaHalfWidthMicrons  & ...
                    coneLocsHexGridY-halfConeApertureMicrons >= yoMicrons-measurementAreaHalfWidthMicrons  & ...
                    coneLocsHexGridY+halfConeApertureMicrons <= yoMicrons+measurementAreaHalfWidthMicrons ));
                if (conesWithin == 0)
                   densityMap(iYpos,iXpos) = nan;
                else
                    densityMap(iYpos,iXpos) = conesWithin / measurementAreaInMM2;
                end
            end
        end
    else
        [X,Y] = meshgrid(gridXPos, gridYPos);
        eccInMeters = sqrt(X.^2 + Y.^2);
        ang = atan2(Y, X)/pi*180;
        [~, ~, densities] = coneSize(eccInMeters(:),ang(:));
        densityMap = reshape(densities, size(X));
    end
    [densityMapSupportX, densityMapSupportY] = meshgrid(gridXPos, gridYPos);
end