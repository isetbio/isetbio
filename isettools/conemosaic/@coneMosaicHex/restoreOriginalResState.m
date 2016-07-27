function restoreOriginalResState(obj)
    
    obj.patternSampleSize = obj.patternSampleSizeOriginatingRectGrid;
    obj.mosaicSize = size(obj.patternOriginatingRectGrid);
    obj.pattern = obj.patternOriginatingRectGrid;    
    obj.setSizeToFOV(obj.fovOriginatingRectGrid);

end