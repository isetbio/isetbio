function cData = contourDataFromDensityMap(spatialSupportXY, zData, zLevels)
    xSupport = squeeze(spatialSupportXY(:,1));
    ySupport = squeeze(spatialSupportXY(:,2));
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;

    iContour = 0;
    levelsNum = 0;
    while (startPoint < dataPoints)
        levelsNum = levelsNum + 1;
        if (levelsNum > 1)
            theLevel = C(1,startPoint);
        end
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        iContour = iContour+1;
        cData{iContour} = struct('faces', f, 'vertices', v);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end