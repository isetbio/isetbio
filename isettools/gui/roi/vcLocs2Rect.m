function rect = vcLocs2Rect(roiLocs)
% Convert roi locs to rect format
%
% Syntax:
%   rect = vcLocs2Rect(roiLocs)
%
% Description:
%    This function takes a Nx2 Matrix entry and returns a 1x4 Vector.
%    rect is a 1x4 spec of a rectangle following the format
%
%       [colMin, rowMin, (cWidth - 1), (rWidth - 1)]
%
%    roiLocs is an Nx2 matrix of (r, c) values
%
%    The code below contains examples of function usage. To access, type
%    'edit vcLocs2Rect.m' into the Command Window.
%
% Inputs:
%    roiLocs - Matrix. Nx2 Matric of (r, c) values.
%
% Outputs:
%    rect    - Vector. 1x4 Vector following the convention described in the
%              description section above.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    ieGetXYCoords, ieRoi2Locs, vcROISelect, vcLineSelect, vcPointSelect
%

% History:
%    xx/xx/12       (c) Imageval Consulting, LLC 2012
%    05/09/18  jnm  Formatting

% Examples:
%{
    roiLocs = vcRect2Locs([2 4 9 7]);
    rect = vcLocs2Rect(roiLocs)
    roiLocs2 = vcRect2Locs(rect);
    isequal(roiLocs, roiLocs2)
%}

if size(roiLocs) ~= 2
    error('Expecting roiLocs as Nx2');
end

rect = zeros(1, 4);

% cMin and rMin
rect(1) = min(roiLocs(:, 2));
rect(2) = min(roiLocs(:, 1));

% cWidth and rWidth
rect(3) = max(roiLocs(:, 2)) - rect(1);
rect(4) = max(roiLocs(:, 1)) - rect(2);

end
