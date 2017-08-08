function squareApertureSize = sizeForSquareApertureFromDiameterForCircularAperture(circularApertureDiameter)
% squareApertureDiameter = sizeForSquareApertureFromDiameterForCircularAperture(circularApertureDiameter)
% 
% Given the diameter of a circular aperature, find the width of a
% square aperture that produces the same area. We use this because in
% isetbio, we specify photoreceptor diameter by height and width of a square
% collecting aperature that are then multiplied to give us the area.
% 
% See also diameterForCircularApertureFromWidthForSquareAperture
%

    area = pi * (circularApertureDiameter/2)^2;
    squareApertureSize = sqrt(area);
end

