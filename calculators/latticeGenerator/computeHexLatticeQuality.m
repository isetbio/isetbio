function [minQValue, qValues] =  computeHexLatticeQuality(rfPositions, triangles)
    
    trianglesNum = size(triangles,1);
    X = rfPositions(:,1);
    Y = rfPositions(:,2);
    
    qValues = zeros(1,trianglesNum);
    for triangleIndex = 1:trianglesNum
        for node = 1:3
            x(node) = X(triangles(triangleIndex,node));
            y(node) = Y(triangles(triangleIndex,node));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        qValues(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end
    
    pointEightPercent = 0.8;
    minQValue = prctile(qValues,pointEightPercent);
end
