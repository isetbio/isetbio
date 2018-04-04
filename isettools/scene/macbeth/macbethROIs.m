function patchLocs = macbethROIs(currentLoc, delta)
% Derive a rectangular region of interest for an MCC chart
%
% Syntax:
%	patchLocs = macbethROIs(currentLoc, [delta])
%
% Description:
%    Find all the locations for a patch centered at the currentLoc and with
%    a spacing of delta. The format of a rectangle is
%       (colMin, rowMin, width, height).
%       The patchLocs is a matrix of N (row, col).
%
% Inputs:
%    currentLoc - The center of the patch.
%    delta      - (Optional) The patch spacing. Default is 10.
%
% Outputs:
%    patchLocs  - A matrix of N (row, col) containing the requested
%                 rectangular region.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    macbethRectangles, macbethDrawRects, macbethSensorValues
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    01/31/18  jnm  Formatting

if notDefined('currentLoc'), error('current location in MCC required'); end
if notDefined('delta'), delta = 10; end  % Get a better algorithm for size

rect(1) = currentLoc(2) - round(delta / 2);
rect(2) = currentLoc(1) - round(delta / 2);
rect(3) = delta;
rect(4) = delta;

patchLocs = ieRoi2Locs(rect);

end