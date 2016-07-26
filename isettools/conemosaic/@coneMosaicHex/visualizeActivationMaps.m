function visualizeActivationMaps(obj, activation)
% Visualize activation maps images for the hex mosaic (all cones +  LMS submosaics)
%
% NPC, ISETBIO TEAM, 2015

    % Compute activation image maps
    [activationImage, activationImageLMScone, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = obj.computeActivationImage(activation);
    activatioRange = [min(activation(:)) max(activation(:))];

    hFig = figure();
    set(hFig, 'Position', [10 10 920 875], 'Color', [1 1 1], 'MenuBar', 'None');
    
    subplotPositions = [...
        0.05 0.53 0.43 0.43; ...
        0.53 0.53 0.43 0.43; ...
        0.05 0.04 0.43 0.43; ...
        0.53 0.04 0.43 0.43];
    
    for subplotIndex = 1:4
        subplot('Position', subplotPositions(subplotIndex,:));
        switch subplotIndex
            case 1
                activationMapImage  = squeeze(activationImageLMScone(:,:,1));
                subplotTitle = 'L-cone submosaic activation';
                showXticks = false;
                showYticks = true;
            case 2
                activationMapImage  = squeeze(activationImageLMScone(:,:,2));
                subplotTitle = 'M-cone submosaic activation';
                showXticks = false;
                showYticks = false;
            case 3
                activationMapImage  = squeeze(activationImageLMScone(:,:,3));
                subplotTitle = 'S-cone submosaic activation';
                showXticks = true;
                showYticks = true;
            case 4
                activationMapImage  = activationImage;
                subplotTitle = 'total mosaic activation';
                showXticks = true;
                showYticks = false;
        end
        
        imagesc(sampledHexMosaicXaxis*1e6, sampledHexMosaicYaxis*1e6, activationMapImage);
        axis 'image'; axis 'xy';
        set(gca, 'CLim', activatioRange);
        if (~showXticks)
            set(gca, 'XTick', []);
        end
        if (~showYticks)
            set(gca, 'YTick', []);
        end
        set(gca, 'FontSize', 14);
        title(subplotTitle, 'FontSize', 16);
    end

    colormap(gray(1024));
    drawnow
end
