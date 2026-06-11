function psfTicksArcMin = spatialMapTicksArcMin(spatialMapRangeArcMin)

    psfTicksMin = (-180:10:180);
   
    if (spatialMapRangeArcMin <= 2)
        psfTicksArcMin  = (-3:0.5:3);
    elseif (spatialMapRangeArcMin <= 5)
        psfTicksArcMin = 0.1*psfTicksMin;
    elseif (spatialMapRangeArcMin <= 10)
        psfTicksArcMin = 0.4*psfTicksMin;
    elseif (spatialMapRangeArcMin <= 20)
        psfTicksArcMin = 1*psfTicksMin;
    elseif (spatialMapRangeArcMin <= 40)
        psfTicksArcMin = 2*psfTicksMin;
    elseif (spatialMapRangeArcMin <= 50)
        psfTicksArcMin = 3*psfTicksMin; 
    elseif (spatialMapRangeArcMin <= 60)
        psfTicksArcMin  = 4*psfTicksMin; 
    elseif (spatialMapRangeArcMin <= 100)
        psfTicksArcMin = 6*psfTicksMin; 
    elseif (spatialMapRangeArcMin<= 200)
        psfTicksArcMin = 10*psfTicksMin;
    elseif (spatialMapRangeArcMin <= 400)
        psfTicksArcMin = 20*psfTicksMin; 
    end
end
