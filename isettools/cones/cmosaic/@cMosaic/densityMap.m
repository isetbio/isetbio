function densityMap2D = densityMap(rfPositions, rfSpacings, sampledPositions)
% Compute a 2D density map from the RF positions
%
% Syntax:
%   [coneDensityMap, coneDensitySupport] = densityMap(obj);
%
% Description:
%    Compute a 2D density map of the cone mosaic 
%
% Inputs:
%    rfPositions          - [Nx2] matrix of (x,y) positions
%    rfSpacings           - [Nx1] matricx of local spacings of the N rfs
%    sampledPositions     - Vectors of sampled positions
%
% Outputs:
%
% Optional key/value pairs:
%    'densityMap2D'     - 2D density maps of cones
%    'densityMapSupport' - spatial support of cone density   

    xSupport = sampledPositions{1};
    ySupport = sampledPositions{2};
    cols = numel(xSupport);
    rows = numel(ySupport);
    densityMap2D = zeros(rows, cols);
    
    for iRow = 1:rows
        yPos = ySupport(iRow);
        [~, idx] = sort(abs(yPos - ySupport));
        yRadius = 0.5*abs(yPos-ySupport(idx(2)));
        for iCol = 1:cols
            xPos = xSupport(iCol);
            [~, idx] = sort(abs(xPos - xSupport));
            xRadius = 0.5*abs(xPos-xSupport(idx(2)));
            nearbyRFindices = find( ...
                (rfPositions(:,1) >= xPos - xRadius) & ...
                (rfPositions(:,1) <= xPos + xRadius) & ...
                (rfPositions(:,2) >= yPos - yRadius) & ...
                (rfPositions(:,2) <= yPos + yRadius));
            meanSpacing = mean(rfSpacings(nearbyRFindices));
            densityMap2D(iRow, iCol) = RGCmodels.Watson.convert.spacingToDensityForHexGrid(meanSpacing);
        end
    end
end

