function restoreOriginalResState(obj)
    
    obj.patternSampleSize = obj.originalResPatternSampleSize;
    obj.mosaicSize = size(obj.originalResPattern);
    obj.pattern = obj.originalResPattern;    
    obj.setSizeToFOV(obj.originalResFOV);

end