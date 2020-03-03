function [conesToMRGCratio, spatialSupport, xLabelString, yLabelString, ratioLabel, ...
            meridianConeToMRGratio, eccUnits] = ...
            compute2DConeToMRGCRFRatio(obj, eccDegsInREVisualSpace,  theReturnedView)
        

   [coneDensity2DMap, coneMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DConeRFDensity(eccDegsInREVisualSpace, theReturnedView);
    
    
   [mRGCDensity2DMap, mRGCMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DmRGCRFDensity(eccDegsInREVisualSpace, theReturnedView);
    
    %On/OFF mRGC Density assuming equal numerosities of ON/OFF
    mRGCDensity2DMap = 0.5*mRGCDensity2DMap;
    mRGCMeridianDensities.temporal = 0.5 * mRGCMeridianDensities.temporal;
    mRGCMeridianDensities.nasal = 0.5 * mRGCMeridianDensities.nasal;
    mRGCMeridianDensities.superior = 0.5 * mRGCMeridianDensities.superior;
    mRGCMeridianDensities.inferior = 0.5 * mRGCMeridianDensities.inferior;
    
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

        