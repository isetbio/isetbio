function [f,v] = facesAndVertices(positions, spacings, shapeOutline)
    shapeSamplesNum = size(shapeOutline,1);
    rfsNum = size(positions, 1);
    X = zeros(rfsNum*shapeSamplesNum,1);
    Y = X;
    f = zeros(rfsNum,shapeSamplesNum);
    for iRF = 1:rfsNum
        ii = (iRF-1)*shapeSamplesNum + (1:shapeSamplesNum);
        f(iRF,:) = ii;
        X(ii,:) = positions(iRF,1) + spacings(iRF) * shapeOutline(:,1);
        Y(ii,:) = positions(iRF,2) + spacings(iRF) * shapeOutline(:,2);
    end
    v = [X(:) Y(:)];
end