function circularApertureDiameter = ...
    diameterForCircularApertureFromWidthForSquareAperture(...
    squareApertureWidth)
% Calculate circular diameter of equal area given square aperture's width
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
%    Examples are contained in the code. To access, type 'edit
%    diameterForCircularApertureFromWidthForSquareAperture.m' in the
%    Command Window.
%
% Input:
%    squareApertureWidth      - Width of square aperture
%
% Output: 
%    circularApertureDiameter - Diameter of circular aperture
%
% Optional key/value pairs:
%    None
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
    area = squareApertureWidth .^ 2;
    circularApertureDiameter = 2.0 * sqrt(area / pi);
end
