function contourCell = plotContours(obj, spacing, col)
% Plots the contours of a spatial receptive field center and surround at a 
% particular standard deviation of the response.
% 
% Inputs: an rgcMosaic object.
% 
% Outputs: contourCell, a cell array containing the RF center coordinates
% in contourCell{:,:,1} and the RF surround coordinates in contourCell{:,:,2}.
% 
% Example:
% 
% % Generate spatial RF contours
% for cellTypeInd = 1:length(obj.mosaic)
%     spatialRFcontours{:,:,:,cellTypeInd} = plotContours(obj.mosaic{cellTypeInd});
% end
% 
% (c) isetbio
% 09/2015 JRG

nCells = size(obj.sRFcenter);

% extent = .5*round(size(obj.sRFcenter{1,1},1)/obj.rfDiameter);
metersPerPixel = (spacing/1e-6)/col;
rfPixels = obj.rfDiameter/metersPerPixel;
extent = .5*round(size(obj.sRFcenter{1,1},1)/rfPixels);

figure;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        % spatialRFFill{xcell,ycell}  = find(obj.sRFcenter{xcell,ycell} > obj.rfDiaMagnitude{xcell,ycell});
        
        hold on;
        [cc,h] = contour(obj.sRFcenter{xcell,ycell},[obj.rfDiaMagnitude{xcell,ycell,1} obj.rfDiaMagnitude{xcell,ycell,1}]);% close;
        cc = bsxfun(@plus,cc,obj.cellLocation{xcell,ycell}' - [1; 1]*extent*rfPixels);
        %         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        contourCell{xcell,ycell,1} = cc;
        
        clear cc h
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr
        [cc,h] = contour(obj.sRFcenter{xcell,ycell},[obj.rfDiaMagnitude{xcell,ycell,2} obj.rfDiaMagnitude{xcell,ycell,2}]);% close;
        %         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
         cc = bsxfun(@plus,cc,obj.cellLocation{xcell,ycell}' - [1; 1]*extent*rfPixels);
        contourCell{xcell,ycell,2} = cc;
        
    end
end
close;