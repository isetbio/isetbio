function contourCell = plotContours(obj)

nCells = size(obj.sRFcenter);

extent = .5*round(size(obj.sRFcenter{1,1},1)/obj.rfDiameter);

figure;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        % spatialRFFill{xcell,ycell}  = find(obj.sRFcenter{xcell,ycell} > obj.rfDiaMagnitude{xcell,ycell});
        
        hold on;
        [cc,h] = contour(obj.sRFcenter{xcell,ycell},[obj.rfDiaMagnitude{xcell,ycell,1} obj.rfDiaMagnitude{xcell,ycell,1}]);% close;
        cc = bsxfun(@plus,cc,obj.cellLocation{xcell,ycell}' - [1; 1]*extent*obj.rfDiameter);
        %         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        contourCell{xcell,ycell,1} = cc;
        
        clear cc h
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr
        [cc,h] = contour(obj.sRFcenter{xcell,ycell},[obj.rfDiaMagnitude{xcell,ycell,2} obj.rfDiaMagnitude{xcell,ycell,2}]);% close;
        %         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
         cc = bsxfun(@plus,cc,obj.cellLocation{xcell,ycell}' - [1; 1]*extent*obj.rfDiameter);
        contourCell{xcell,ycell,2} = cc;
        
    end
end
close;