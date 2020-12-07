function sizeMicrons = sizeVisualDegsToSizeRetinalMicrons(sizeVisualDegs, eccDegs)
% Convert a retinal size in visual degrees at given ecc (also in degs) to its
% equivalent linear size in retinal microns
%
% Syntax:
%   sizeMicrons = sizeVisualDegsToSizeRetinalMicrons(sizeVisualDegs, eccDegs)
%
    
    if (numel(eccDegs) == 2) && (numel(sizeVisualDegs) > 1)
        % Treat eccDegs as a vector of (x,y) eccs
        eccDegs = sqrt(sum(eccDegs.^2,2));
    elseif (numel(eccDegs) == 1) || (numel(sizeVisualDegs) == 1)
        % Do nothing
    else
        error('dont know how to treat the passed dimensions');
    end
    
    leftCoordDegs  = (eccDegs - sizeVisualDegs/2);
    rightCoordDegs = (eccDegs + sizeVisualDegs/2);
    sizeMicrons =  1e3 * ( RGCmodels.Watson.convert.rhoDegsToMMs(rightCoordDegs) ...
                          -RGCmodels.Watson.convert.rhoDegsToMMs(leftCoordDegs));
          
end

