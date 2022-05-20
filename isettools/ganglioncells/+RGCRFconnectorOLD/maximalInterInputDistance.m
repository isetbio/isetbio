function d = maximalInterInputDistance(coneRFpos)
    d = pdist2(coneRFpos, coneRFpos);
    d = max(d(:));
end