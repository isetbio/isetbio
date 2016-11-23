function hex2Dmap = reshapeHex3DmapToHex2Dmap(obj, hex3Dmap)

nonNullCones = obj.pattern(obj.pattern>1);

iLdest = find(nonNullCones==2);
iMdest = find(nonNullCones==3);
iSdest = find(nonNullCones==4);

iLsource = find(obj.pattern==2);
iMsource = find(obj.pattern==3);
iSsource = find(obj.pattern==4);

timePointsNum = size(hex3Dmap,3);
hex3Dmap = reshape(hex3Dmap, [size(obj.pattern,1)*size(obj.pattern,2) timePointsNum]);
hex2Dmap = zeros(numel(nonNullCones), timePointsNum);
hex2Dmap(iLdest,:) = hex3Dmap(iLsource,:);
hex2Dmap(iMdest,:) = hex3Dmap(iMsource,:);
hex2Dmap(iSdest,:) = hex3Dmap(iSsource,:);

end

