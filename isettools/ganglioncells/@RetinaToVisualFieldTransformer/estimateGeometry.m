function [theCentroid, theAxesLengths, theRotationAngle] = estimateGeometry(supportX, supportY, zData)
    % Compute orientation, centroid, and major/minor axis lengths
    binaryImage = zData;
    m1 = min(binaryImage(:));
    m2 = max(binaryImage(:));
    binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
    s = regionprops('table', binaryImage,'Orientation', 'Centroid', 'MinorAxisLength', 'MajorAxisLength');
    theCentroid = s.Centroid(1,:);
    theMinorAxisLength = s.MinorAxisLength(1);
    theMajorAxisLength = s.MajorAxisLength(1);
    theRotationAngle = s.Orientation(1);

    % The computed centroid and axis lengths are in pixels. Convert them to units of spatial support
    % to serve as initial Gaussian parameter values
    xx = [round(theCentroid(1)-theMajorAxisLength*0.5) round(theCentroid(1)+theMajorAxisLength*0.5)];
    yy = [round(theCentroid(2)-theMinorAxisLength*0.5) round(theCentroid(2)+theMinorAxisLength*0.5)];

    xx(1) = max([1 xx(1)]);
    xx(2) = min([numel(supportX) xx(2)]);

    yy(1) = max([1 yy(1)]);
    yy(2) = min([numel(supportY) yy(2)]);

    theAxesLengths(1) = (supportY(yy(2)) - supportY(yy(1)));
    theAxesLengths(2) = (supportX(xx(2)) - supportX(xx(1)));
    theCentroid(1) = supportX(round(theCentroid(1)));
    theCentroid(2) = supportY(round(theCentroid(2)));
end