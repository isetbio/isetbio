function sizeMicrons = sizeVisualDegsToSizeRetinalMicrons(obj,sizeVisualDegs, eccDegs)
% Convert a retinal size in visual degrees at given ecc (also in degs) to its
% equivalent linear size in retinal microns
%
    
    if (numel(eccDegs) == 2) && (numel(sizeVisualDegs) > 1)
        % Treat eccDegs as a vector of (x,y) eccs
        eccDegs = sqrt(sum(eccDegs.^2,2));
    elseif (numel(eccDegs) == 1) || (numel(sizeVisualDegs) == 1)
        % Do nothing
    elseif (size(sizeVisualDegs,1) == size(eccDegs,1)) &&  (size(sizeVisualDegs,2)*size(eccDegs,2) == 1)
        % Do nothing
    elseif (size(sizeVisualDegs,2) == size(eccDegs,2)) &&  (size(sizeVisualDegs,1)*size(eccDegs,1) == 1)
        % Do nothing
    else
        error('dont know how to treat the passed dimensions');
    end
    
    leftCoordDegs  = (eccDegs - sizeVisualDegs/2);
    rightCoordDegs = (eccDegs + sizeVisualDegs/2);
    sizeMicrons =  1e3 * ( obj.customDegsToMMsConversionFunction(rightCoordDegs) ...
                          -obj.customDegsToMMsConversionFunction(leftCoordDegs));
          
end
