function [retinalConeActivations, visualConeActivations, ...
    withinConeAperturesRetinalConeImage, ...
    withinConeAperturesVisualConeImage, withinConeAperturesGrid] = integrateConeImageWithinConeAperturesForDeconvolution(...
    retinalConeImage, visualConeImage, coneMask, conePosDegs, coneAperturesDegs, thePSFsupportDegs)

    conesNum = size(conePosDegs,1);
    retinalConeActivations = zeros(1, conesNum);
    visualConeActivations = zeros(1, conesNum);
    [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
    
    retinalConeImage = retinalConeImage .* coneMask;
    visualConeImage = visualConeImage .* coneMask;
    
    allPixelsWithinConeApertures = [];
    
    coneOutline.x = 0.5*cosd(0:2:360);
    coneOutline.y = 0.5*sind(0:2:360);
    
    for iCone = 1:conesNum
        xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*coneOutline.x;
        yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*coneOutline.y;
        [in,on] = inpolygon(X(:), Y(:), xr, yr);
        withinAperturePixelIndices = find((in==true));
        allPixelsWithinConeApertures = cat(1,allPixelsWithinConeApertures, withinAperturePixelIndices(:));
        if (numel(withinAperturePixelIndices) == 0)
            retinalConeActivations(iCone) = 0;
            visualConeActivations(iCone) = 0;
        else
            netRetinalConeActivation = sum(retinalConeImage(withinAperturePixelIndices));
            netVisualConeActivation  = sum(visualConeImage(withinAperturePixelIndices));
            retinalConeActivations(iCone) = netRetinalConeActivation;
            visualConeActivations(iCone) = netVisualConeActivation;
        end
    end

    
    withinConeAperturesRetinalConeImage = retinalConeImage(:);
    withinConeAperturesRetinalConeImage = withinConeAperturesRetinalConeImage(allPixelsWithinConeApertures);
    
    withinConeAperturesVisualConeImage = visualConeImage(:);
    withinConeAperturesVisualConeImage = withinConeAperturesVisualConeImage (allPixelsWithinConeApertures);
    
    withinConeAperturesGrid = [X(:) Y(:)];
    withinConeAperturesGrid = withinConeAperturesGrid(allPixelsWithinConeApertures,:);
    
end