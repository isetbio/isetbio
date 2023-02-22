function RGCindices = rgcIndicesAlongMeridianWithAngle(...
    theMeridianAngle, theMeridianRadius, rgcRFpositionsDegs, maxRGCsNum)
        
    x = theMeridianRadius * cosd(theMeridianAngle);
    y = theMeridianRadius * sind(theMeridianAngle);

    theROI = regionOfInterest('shape', 'line', 'from', [x,y], 'to', [-x,-y], 'thickness', 0.01);
    samplingPoints = 400;  % sample the perimeter of the ROI along 1000 points');
    pointsPerSample = 10;  % return up to 10 points for each sample along the perimeter');
    maxDistance = 0.1;     % points must be no further than 0.1 degs away from the closest perimeter sample');
    idx = theROI.indicesOfPointsAround(...
        rgcRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
    
    if (~isempty(maxRGCsNum))
        skip = floor(numel(idx)/maxRGCsNum);
    else
        skip = 1;
    end
    idx = idx(1:skip:numel(idx));
    RGCindices = idx;

    % Sort the indices
    if (theMeridianAngle == 90) || (theMeridianAngle == 270)
        % Sort according to y-coord
        [~,idx] = sort(squeeze(rgcRFpositionsDegs(RGCindices,2)), 'ascend');
    else
        % Sort according to x-coord
        [~,idx] = sort(squeeze(rgcRFpositionsDegs(RGCindices,1)), 'ascend');
    end
    RGCindices = RGCindices(idx);

end
