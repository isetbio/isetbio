function sizeDegs = sizeRetinalMicronsToSizeVisualDegs(obj, sizeRetinalMicrons, eccMicrons)
% Convert a retinal size in microns at given ecc (also in microns) to its
% equivalent angular size in visual degrees
    
    if (numel(eccMicrons) == 2) && (numel(sizeRetinalMicrons) > 1)
        % Treat eccMicrons as a vector of (x,y) eccs
        eccMicrons = sqrt(sum(eccMicrons.^2,2));
    elseif (numel(eccMicrons) == 1) || (numel(sizeRetinalMicrons) == 1)
         % Do nothing
    elseif (size(sizeRetinalMicrons,1) == size(eccMicrons,1)) && (size(sizeRetinalMicrons,2)*size(eccMicrons,2) == 1)
        % Do nothing
    elseif (size(sizeRetinalMicrons,2) == size(eccMicrons,2)) && (size(sizeRetinalMicrons,1)*size(eccMicrons,1) == 1)
        % Do nothing
    else
         error('dont know how to treat the passed dimensions');
    end
    
    leftCoordMM  = (eccMicrons - sizeRetinalMicrons/2)/1e3;
    rightCoordMM = (eccMicrons + sizeRetinalMicrons/2)/1e3;
    sizeDegs =  obj.customMMsToDegsConversionFunction(rightCoordMM) ...
               -obj.customMMsToDegsConversionFunction(leftCoordMM);
          
end

