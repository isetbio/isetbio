function circularApertureDiameter = ...
    diameterForCircularApertureFromWidthForSquareAperture(...
    squareApertureWidth)
% Calculate diameter of circular aperture of identical area given width for
% square aperture
%
% Syntax:
%    circularApertureDiameter = ...
%       diameterForCircularApertureFromWidthForSquareAperture(...
%       squareApertureWidth)
% 
% Description
%    Given the width of a square aperature, find the diameter of a circular
%    aperture that produces the same area. We use this because in isetbio,
%    we specify photoreceptor diameter by height and width of a square
%    collecting aperature that are then multiplied to give us the area.
%
% Input:
%    squareApertureWidth      - Width of square aperture
%
% Output: 
%    circularApertureDiameter - Diameter of circular aperture
%
% See Also:
%    sizeForSquareApertureFromDiameterForCircularAperture
%

% Example:
%{
    circularApertureDiameter = ...
        diameterForCircularApertureFromWidthForSquareAperture(...
        squareApertureWidth)
%}
    area = squareApertureWidth^2;
    circularApertureDiameter = 2.0*sqrt(area/pi);
end
