function radius = selectPlotSupport(data, prct, minRadius)
% Used with getMiddleMatrix to pull out the 'interesting' center of a plot
%
% Syntax:
%   sz = selectPlotSupport(data, [prct = 0.01],[minRadius = 25])
%
% Description:
%   Often we have a large surface to plot but the interesting
%   (non-zero) part is near the middle of the matrix. Rather than
%   plotting the entire surf or mesh(data) we pull out a central
%   region so it is easier to see the relevant data.
%
%   This the method chooses the radius of a central region by finding
%   values that exceed 1 percent (adjusted by the prct variable) of
%   the max value in the central row.  A minimum radius is 25 samples,
%   also adjustable by a parameter.
%
%   The method is used in conjunction with getMiddleMatrix.
%
% Inputs:
%    data - The data set to plot
%    prct - (Optional) What percentage of the middle you wish to display,
%           in decimal form. Default is 0.01 (1%)
%
% Outputs:
%    sz   - Size (radius) of the extracted region
%
% Notes:
%	 * [Note: XXX - What if data are a vector? Can we adjust this routine
%	   to make it work?]
%
%  See Also:
%    oiPlot, getMiddleMatrix

% Examples:
%{
 g = fspecial('gaussian',256,30);
 g = ieScale(g,0,1);
 r = selectPlotSupport(g,0.1);
 vcNewGraphWin; imagesc(g); axis image
 ieShape('circle','center',[128,128],'radius',r,'color','white');
 colorbar;
%}

%%

if notDefined('prct'), prct = 0.01; end
if notDefined('minRadius'), minRadius = 25; end

r  = size(data, 1);
mx = max(data(:));
centerRow = round(r / 2);

%% Find the location in the center row at the edge of the prct max
l = (data(centerRow, :) < prct * mx);

if max(l) == 0, radius = centerRow - 1;
else
    % The returned value is at least 25.
    [~, idx] = max(data(centerRow, l));
    radius = centerRow - idx;
    if radius < 25
        warning('Selected radius is very small. Increasing to %d',minRadius);
    end
end

end