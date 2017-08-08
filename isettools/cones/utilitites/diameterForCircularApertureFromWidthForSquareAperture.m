function circularApertureDiameter = diameterForCircularApertureFromWidthForSquareAperture(squareApertureWidth)
% circularApertureDiameter = diameterForCircularApertureFromWidthForSquareAperture(squareApertureWidth)
% 
% Given the width of a square aperature, find the diameter of a
% circular aperture that produces the same area. We use this because in
% isetbio, we specify photoreceptor diameter by height and width of a square
% collecting aperature that are then multiplied to give us the area.
%
% See also sizeForSquareApertureFromDiameterForCircularAperture
%

    area = squareApertureWidth^2;
    circularApertureDiameter = 2.0*sqrt(area/pi);
end
