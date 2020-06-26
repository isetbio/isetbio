function contourCell = plotContours(obj, spacing, col)
% Plots spatial receptive field center and surround contours
%
% Syntax:
%   contourCell = plotContours(obj, spacing, col)
%
% Description:
%    Plot a spatial receptive field with center and surround contours.
%
%    This function contains inline examples of usage. In order to access
%    these, type 'edit plotContours.m' into the Command Window.
%
% Inputs:
%   obj         - Object. An rgcMosaic object.
%   spacing     - Apparently, patch size of the inner retina is meant
%   col         - Number of columns in the mosaic
%
% Outputs:
%   contourCell - Cell. A cell array containing the RF coordinates, with
%                 contourCell{:, :, 1} denoting the center, and 
%                 contourCell{:, :, 2} denoting the surround.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/15  JRG  (c) isetbio team, 2015
%    06/06/19  JNM  Documentation pass

% Example:
%{
    for cellTypeInd = 1:length(obj.mosaic)
        spatialRFcontours{:, :, :, cellTypeInd} = ...
            plotContours(obj.mosaic{cellTypeInd});
    end
%}

nCells = size(obj.sRFcenter);

% extent = .5 * round(size(obj.sRFcenter{1, 1}, 1) / obj.rfDiameter);
metersPerPixel = (spacing / 1e-6) / col;
rfPixels = obj.rfDiameter / metersPerPixel;
extent = .5 * round(size(obj.sRFcenter{1, 1}, 1) / rfPixels);
vcNewGraphWin;
% Replace this code with multiple circle draws, I think.
contourCell = cell(nCells(1), nCells(2), 2);
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        hold on;
        cc = contour(obj.sRFcenter{xcell, ycell}, ...
            [obj.rfDiaMagnitude{xcell, ycell, 1}, ...
            obj.rfDiaMagnitude{xcell, ycell, 1}]);  % close;
        cc = bsxfun(@plus, cc, ...
            obj.cellLocation{xcell, ycell}' - [1; 1] * extent * rfPixels);
        cc(:, 1) = [NaN; NaN];
        contourCell{xcell, ycell, 1} = cc;

        clear cc h

        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr
        cc = contour(obj.sRFcenter{xcell, ycell}, ...
            [obj.rfDiaMagnitude{xcell, ycell, 2}, ...
            obj.rfDiaMagnitude{xcell, ycell, 2}]);  % close;
        cc(:, 1) = [NaN; NaN];
        cc = bsxfun(@plus, cc, ...
            obj.cellLocation{xcell, ycell}' - [1; 1] * extent * rfPixels);
        contourCell{xcell, ycell, 2} = cc;
    end
end
close;

end
