function [f,v] = facesAndVertices(positions, spacings, diskOutline)
    thetaSamples = size(diskOutline,1);
    rfsNum = size(positions, 1);
    X = zeros(rfsNum*thetaSamples,1);
    Y = X;
    f = zeros(rfsNum,thetaSamples);
    for iRF = 1:rfsNum
        ii = (iRF-1)*thetaSamples + (1:thetaSamples);
        f(iRF,:) = ii;
        X(ii,:) = positions(iRF,1) + spacings(iRF) * diskOutline(:,1);
        Y(ii,:) = positions(iRF,2) + spacings(iRF) * diskOutline(:,2);
    end
    v = [X(:) Y(:)];
end
