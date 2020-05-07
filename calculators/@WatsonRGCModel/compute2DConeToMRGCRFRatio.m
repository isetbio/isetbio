function [conesToMRGCratio, spatialSupport, xLabelString, yLabelString, ratioLabel, ...
            meridianConeToMRGratio, eccUnits] = ...
            compute2DConeToMRGCRFRatio(obj, eccDegsInREVisualSpace,  theReturnedView)
        

   [coneDensity2DMap, coneMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DConeRFDensity(eccDegsInREVisualSpace, theReturnedView);
    
    
   [mRGCDensity2DMap, mRGCMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DmRGCRFDensity(eccDegsInREVisualSpace, theReturnedView, 'subtype', 'ON');
    
    conesToMRGCratio = coneDensity2DMap./mRGCDensity2DMap;
    meridianConeToMRGratio.ecc = coneMeridianDensities.ecc;
    meridianConeToMRGratio.temporal = coneMeridianDensities.temporal ./ mRGCMeridianDensities.temporal;
    meridianConeToMRGratio.nasal    = coneMeridianDensities.nasal ./ mRGCMeridianDensities.nasal;
    meridianConeToMRGratio.superior = coneMeridianDensities.superior ./ mRGCMeridianDensities.superior;
    meridianConeToMRGratio.inferior = coneMeridianDensities.inferior ./ mRGCMeridianDensities.inferior;
    
    ratioLabel = 'cones to mRGC ratio';
    spatialSupport = densitySupport;
    xLabelString = horizontalMeridianLabel;
    yLabelString = verticalMeridianLabel;
    eccUnits = supportUnits;
end

        