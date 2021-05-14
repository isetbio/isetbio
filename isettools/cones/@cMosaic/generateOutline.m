% Static method to generate an ROI outline from an ROIstruct
function roiOutline = generateOutline(roi)
    if (strcmp(roi.shape, 'rect'))
        roiOutline.x = roi.center(1) + 0.5*roi.width  * [-1 -1  1  1 -1];
        roiOutline.y = roi.center(2) + 0.5*roi.height * [-1  1  1 -1 -1];
    else
        a = 0.5*roi.minorAxisDiameter;
        b = 0.5*roi.majorAxisDiameter;
        cosOutline = cosd(0:10:360);
        sinOutline = sind(0:10:360);
        roiOutline.x = roi.center(1) + a*cosOutline*cosd(roi.rotation) - b*sinOutline*sind(roi.rotation);
        roiOutline.y = roi.center(2) + a*cosOutline*sind(roi.rotation) + b*sinOutline*cosd(roi.rotation);
    end
end
