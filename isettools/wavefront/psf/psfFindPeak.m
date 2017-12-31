function [peakRow, peakCol] = psfFindPeak(input)
% Find row/col corresponding to max of two dimensional PSF input.
%
% Syntax:
%   [peakRow, peakCol] = psfFindPeak(input)
%
% Description:
%    Find row/col corresponding to max of two dimensional PSF input.
%
%    This method works reasonably when there is more than one maximum in
%    the input, in that it returns the coordinates of one of the maxima.
%    No thought is investing in figuring out which one if there are
%    multiple locations with the identical maximal value -- it's just
%    whatever the max() routine picks.
%
%    The old method of finding row and column  maxima can screw up in this
%    case, by returing the coordinates of a point that isn't a maximum.
%
% Inputs:
%    input   - Two-dimensional PSF input
%
% Outputs:
%    peakRow - Corresponding row to PSF Max
%    peakCol - Corresponding column to PSF Max
%

% History:
%    12/22/09  dhb  Encapsulate this as a function with improved method.
%    11/10/17  jnm  Comments & formatting

% Examples:
%{
    [pR pC] = psfFindPeak([0 1 10 1; 1 2 5 1; 1 1 1 4])
%}

[m, n] = size(input);
[~, ind] = max(input(:));
[peakRow, peakCol] = ind2sub([m, n], ind);
if (input(peakRow, peakCol) ~= max(input(:)))
    error('Value at max location is not input maximum. Hmmm.');
end

end
