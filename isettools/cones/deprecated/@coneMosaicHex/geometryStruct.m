function cmStruct = geometryStruct(obj)
% % Return a struct with the mosaic geometry
%
% Syntax:
%    cmStruct = geometryStruct(obj)
%
% Description:
%    Return a struct with the mosaic geometry (cone positions, Delaunay 
%    triangles, and cone aperture sizes)
%       
%
% Inputs:
%    obj - The cone mosaic hex object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% History:
%    3/3/19  NPC  ISETBIO TEAM, 2019

    % cone positions in degrees
    cmStruct.coneLocs = obj.coneLocsHexGrid * 1e6 / obj.micronsPerDegree;
    
    % cone types (K/L/M/S)
    cmStruct.coneTypes = obj.coneTypesHexGrid;
    
    % cone aperture in degrees
    if (obj.eccBasedConeDensity)
        % ecc-varying aperture
        cmStruct.coneApertures = (obj.computeApertureDiametersHexGrid()/obj.micronsPerDegree)';
    else
        % non ecc-varying aperture
        cmStruct.coneApertures = obj.customInnerSegmentDiameter * (ones(size(obj.computeApertureDiametersHexGrid()))/obj.micronsPerDegree)';
    end
    
    % Delaunayn triangularization
    xHex = cmStruct.coneLocs(:,1);
    yHex = cmStruct.coneLocs(:,2);
    cmStruct.triangles = delaunayn([xHex(:), yHex(:)]);
    
    plotData = false;
    if (plotData)
    trianglesNum = size(cmStruct.triangles, 1);
    
    figure();
    axis 'equal'
    hold on
    for triangleIndex = 1:trianglesNum 
        coneIndices = cmStruct.triangles(triangleIndex, :);
        xCoords(1:3) = xHex(coneIndices); xCoords(4) = xCoords(1);
        yCoords(1:3) = yHex(coneIndices); yCoords(4) = yCoords(1);
        plot(xCoords, yCoords,'k-');
        if (mod(triangleIndex,1000) == 0)
        drawnow;
        end
    end
    
    % Compute hex-index for all triangles
    hexIndices = computeHexIndex(cmStruct.coneLocs, cmStruct.triangles);
    hexIndicesBins = 0.0:0.01:1.0;
    [counts,centers] = hist(hexIndices, hexIndicesBins);
    
    
    figure()
    bar(centers,counts,1)
    end
    
end

function hexIndices = computeHexIndex(coneLocs, triangles)
    
    trianglesNum = size(triangles,1);
    Xlocs = coneLocs(:,1);
    Ylocs = coneLocs(:,2);
    
    hexIndices = zeros(1,trianglesNum);
    x = zeros(1,3); y = x;
    for triangleIndex = 1:trianglesNum
        for node = 1:3
            x(node) = Xlocs(triangles(triangleIndex,node));
            y(node) = Ylocs(triangles(triangleIndex,node));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        hexIndices(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end

end
