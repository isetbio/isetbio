function XYZ = xyy2xyz(xyy)
% Convert data from CIE xyY to CIE XYZ values
%
% Syntax:
%   XYZ = xyy2xyz(xyy)
%
% Description:
%    It is common to represent the color of a light using xyY (chromaticity
%    coordinates and luminance). This routine converts from the xyY
%    representation to the standard XYZ representation.
%
%    The values are all represented in n X 3 format. For XYZ, the first
%    column is X, second  column is Y and third is Z. Similarly, the input
%    columns must be x,y, and Y.
%
%    Formula:
%       X = (x / y) * Y,
%       Z = ((1 - x - y)/y) * Y
%       Also, note that Y / y = X + Y + Z
%
%    This function contains examples of usage inline. To access these, type
%    'edit xyy2xyz.m' into the Command Window.
%
% Inputs:
%    xyy - Matrix. The chromacity Coordinates and Luminance.
%
% Outputs:
%    XYZ - Matrix. The standard representation of light and color.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    11/01/17  jnm  Comments & formatting
%    11/17/17  jnm  Formatting
%    07/16/19  JNM  Formatting update

% Examples:
%{
    % Use PTB XYZToxyY, included in isetbio, to produce test input
    testXYZs = [[1 2 1]' [2 1 0.5]' [1 1 1]' [0.6 2.3 4]']
    ptbxyYs = XYZToxyY(testXYZs);
    isetXYZs = xyy2xyz(ptbxyYs')'
%}

% Check input
if size(xyy, 2) ~= 3, error('Input must be x,y,Y in the rows.'); end

XYZ = zeros(size(xyy));

% = Y
XYZ(:, 2) = xyy(:, 3);

% X + Y + Z = Y/y
sXYZ = xyy(:, 3) ./ xyy(:, 2);

% X = (x/y)*Y
XYZ(:, 1) = (xyy(:, 1) ./ xyy(:, 2)) .* xyy(:, 3);

% Z = (X + Y + Z) - Y - X
XYZ(:, 3) = sXYZ - xyy(:, 3) - XYZ(:, 1);

end
