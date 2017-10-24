function squareApertureSize = ...
    sizeForSquareApertureFromDiameterForCircularAperture(...
    circularApertureDiameter)
% calculate width of square aperture of identical area given diameter for
% circular aperture
% 
% Syntax:
%       squareApertureWidth = ...
%           sizeForSquareApertureFromDiameterForCircularAperture(...
%           circularApertureDiameter)
% 
% Description:
%    Given the diameter of a circular aperature, find the width of a square
%    aperture that produces the same area. We use this because in isetbio,
%    we specify photoreceptor diameter by height and width of a square
%    collecting aperature that are then multiplied to give us the area.
% 
% Input:
%    circularApertureDiameter - Diameter of circular aperture
%
% Output: 
%    squareApertureSize       - Width of square aperture
%
% See Also:
%    diameterForCircularApertureFromWidthForSquareAperture
%

% Example:
%{
    squareApertureWidth = ...
        sizeForSquareApertureFromDiameterForCircularAperture(...
        circularApertureDiameter)
%}
    area = pi * (circularApertureDiameter/2)^2;
    squareApertureSize = sqrt(area);
end

