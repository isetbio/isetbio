function q = meshQuality(rfLocs, triangles)
    
    trianglesNum = size(triangles,1);
    X = rfLocs(:,1);
    Y = rfLocs(:,2);
    
    q = zeros(1,trianglesNum);
    x = zeros(1,3);
    y = zeros(1,3);
    
    for triangleIndex = 1:trianglesNum
        for vertex = 1:3
            x(vertex) = X(triangles(triangleIndex,vertex));
            y(vertex) = Y(triangles(triangleIndex,vertex));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        q(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end
end


