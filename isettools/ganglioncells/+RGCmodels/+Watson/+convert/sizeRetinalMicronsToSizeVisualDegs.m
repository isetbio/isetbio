function sizeDegs = sizeRetinalMicronsToSizeVisualDegs(sizeRetinalMicrons, eccMicrons)
% Convert a retinal size in microns at given ecc (also in microns) to its
% equivalent angular size in visual degrees
    
    leftCoordMM  = (eccMicrons - sizeRetinalMicrons/2)/1e3;
    rightCoordMM = (eccMicrons + sizeRetinalMicrons/2)/1e3;
    sizeDegs =  RGCmodels.Watson.convert.rhoMMsToDegs(rightCoordMM) ...
               -RGCmodels.Watson.convert.rhoMMsToDegs(leftCoordMM);
          
end

