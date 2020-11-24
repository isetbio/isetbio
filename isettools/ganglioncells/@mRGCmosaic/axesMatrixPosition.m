function axPos = axesMatrixPosition(obj)
    % Compute normalized RGC positions
    mins = min(obj.rgcRFpositionsDegs,[],1);
    maxs = max(obj.rgcRFpositionsDegs,[],1);
    normalizedPos = obj.rgcRFpositionsDegs*0;
    normalizedPos(:,1) = 0.02+0.9*(obj.rgcRFpositionsDegs(:,1)-mins(1))/(maxs(1)-mins(1));
    normalizedPos(:,2) = 0.02+0.9*(obj.rgcRFpositionsDegs(:,2)-mins(2))/(maxs(2)-mins(2));
    
    rgcsNum = size(obj.rgcRFpositionsDegs,1);
    w = 0.6 * mean(obj.rgcRFspacingsDegs) / (maxs(1)-mins(1));
    h = w;
    axPos = zeros(rgcsNum, 4);
    for iRGC = 1:rgcsNum
        axPos(iRGC,:) = [normalizedPos(iRGC,1), ...
                 normalizedPos(iRGC,2), ...
                 w, h];  
    end
end

