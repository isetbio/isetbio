function d = maximalInterInputDistance(RFpos)
    d = MosaicConnector.pdist2(RFpos, RFpos);
    d = max(d(:));
end