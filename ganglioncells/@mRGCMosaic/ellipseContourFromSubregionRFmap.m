function contourData = ellipseContourFromSubregionRFmap(xSupport, ySupport, RF, zLevel, centerSubregionContourSamples)

    % Binarize
    RF = RF / max(RF(:));
    RF(RF<zLevel) = 0.0;
    RF(RF>0) = 1.0;
    BW = imbinarize(RF);

    % Extract the maximum area
    BW = imclearborder(BW);
    BW = bwareafilt(BW,1);

    % Calculate centroid, orientation and major/minor axis length of the ellipse
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    if (isempty(s))
        fprintf(2, 'Could not fit an ellipse to subregion map. Returning empty contour.\n')
        contourData = [];
        return;
    end

    % Calculate the ellipse line
    theta = linspace(0, 2*pi, centerSubregionContourSamples);
    col = (s.MajorAxisLength/2)*cos(theta);
    row = (s.MinorAxisLength/2)*sin(theta);
    M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];

    x = D(1,:);
    y = D(2,:);
    x = (x-1)/(numel(xSupport)-1) * (xSupport(end)-xSupport(1)) + xSupport(1); 
    y = (y-1)/(numel(ySupport)-1) * (ySupport(end)-ySupport(1)) + ySupport(1); 

    v = [x(:) y(:)];
    f = 1:numel(x);
    s = struct('faces', f, 'vertices', v);
    contourData = s;

end

