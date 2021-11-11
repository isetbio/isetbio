function roiOutline = generateOutline(roi)
    if (strcmp(roi.shape, 'rect'))
        a = 0.5;
        b = 0.5;
        cosOutline = [-1 -1  1  1 -1]*roi.width;
        sinOutline = [-1  1  1 -1 -1]*roi.height;
    else
        a = 0.5*roi.minorAxisDiameter;
        b = 0.5*roi.majorAxisDiameter;
        cosOutline = cosd(0:10:360);
        sinOutline = sind(0:10:360);
    end
    roiOutline.x = roi.center(1) + a*cosOutline*cosd(roi.rotation) - b*sinOutline*sind(roi.rotation);
    roiOutline.y = roi.center(2) + a*cosOutline*sind(roi.rotation) + b*sinOutline*cosd(roi.rotation);
end