% Static method to convert an ROIoutline in degs to an ROIoutline in microns
function roiOutlineMicrons = convertOutlineToMicrons(obj, roiOutlineDegs)
    roiOutlineMicrons.x = obj.distanceDegreesToDistanceMicronsForCmosaic(roiOutlineDegs.x);
    roiOutlineMicrons.y = obj.distanceDegreesToDistanceMicronsForCmosaic(roiOutlineDegs.y);
end