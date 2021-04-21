% Static method to convert an ROIoutline in degs to an ROIoutline in microns
function roiOutlineMicrons = convertOutlineToMicrons(roiOutlineDegs,micronsPerDegreeApproximation)
    % Convert degs to microns
    if (~isempty(micronsPerDegreeApproximation))
        % Convert outline from degs to microns using the passed microns/deg approximation
        roiOutlineMicrons.x = roiOutlineDegs.x * micronsPerDegreeApproximation;
        roiOutlineMicrons.y = roiOutlineDegs.y * micronsPerDegreeApproximation;
    else
        roiOutlineMicrons.x = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(roiOutlineDegs.x); 
        roiOutlineMicrons.y = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(roiOutlineDegs.y); 
    end
end