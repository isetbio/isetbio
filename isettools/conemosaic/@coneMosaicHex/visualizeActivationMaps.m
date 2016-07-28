function visualizeActivationMaps(obj, activation, varargin)
% Visualize activation maps images for the hex mosaic (all cones +  LMS submosaics)
%
% NPC, ISETBIO TEAM, 2015

    p = inputParser;
    p.addParameter('figureSize', [920 875], @isnumeric);
    p.addParameter('signalName', ' ', @ischar);
    p.parse(varargin{:});        
            
    % Compute activation image maps
    [activationImage, activationImageLMScone, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = obj.computeActivationImage(activation);
    activeConesActivations = activation(obj.pattern>1);
    %activationRange = prctile(activeConesActivations, [10 90]);
    activationRange = [0 max(activeConesActivations(:))];
    
    hFig = figure();
    set(hFig, 'Position', cat(2, [10 10], p.Results.figureSize), 'Color', [1 1 1]); % , 'MenuBar', 'None');
    
    subplotPositions = [...
        0.02 0.52 0.45 0.45; ...
        0.49 0.52 0.45 0.45; ...
        0.02 0.03 0.45 0.45; ...
        0.49 0.03 0.45 0.45];
    
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
        set(gca, 'CLim', activationRange);
        if (~showXticks)
            set(gca, 'XTick', []);
        end
        if (~showYticks)
            set(gca, 'YTick', []);
        end
        set(gca, 'FontSize', 14);
        title(subplotTitle, 'FontSize', 16);
        
        if (subplotIndex == 4)
            % Add colorbar
            originalPosition = get(gca, 'position');
            hCbar = colorbar('eastoutside', 'peer', gca); % , 'Ticks', cbarStruct.ticks, 'TickLabels', cbarStruct.tickLabels);
            hCbar.Orientation = 'vertical';
            hCbar.Label.String = p.Results.signalName;
            hCbar.FontSize = 14;
            hCbar.FontName = 'Menlo';
            hCbar.Color = [0.2 0.2 0.2];
            % The addition changes the figure size, so undo this change
            newPosition = get(gca, 'position');
            set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
        end
        
    end

    cmap = jet(1024);
    cmap(1,:) = 0;
    colormap(cmap);
    drawnow
end
