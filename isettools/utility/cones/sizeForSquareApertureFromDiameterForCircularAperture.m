function squareApertureSize = ...
    sizeForSquareApertureFromDiameterForCircularAperture(...
    circularApertureDiameter)
% calculate square's width given the diameter of an equal circular aperture
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
%    Examples are contained in the code. To access, type 'edit
%    sizeForSquareApertureFromDiameterForCircularAperture.m' into the
%    Command Window.
%
% Input:
%    circularApertureDiameter - Diameter of circular aperture
%
% Output: 
%    squareApertureSize       - Width of square aperture
%
% Optional key/value pairs:
%    None.
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
    area = pi * (circularApertureDiameter / 2) .^ 2;
    squareApertureSize = sqrt(area);
end