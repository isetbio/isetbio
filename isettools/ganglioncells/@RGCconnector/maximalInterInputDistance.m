function d = maximalInterInputDistance(coneRFpos)
    d = RGCconnector.pdist2(coneRFpos, coneRFpos);
    d = max(d(:));
end