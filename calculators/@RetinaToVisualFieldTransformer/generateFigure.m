function generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
    visualRFparams, retinalRFparams, achievedVisualRFparams, ...
    eccDegs, testSubjectID, maxSpatialSupportDegs, figNo)

    maxPSF = max([max(thePSFData.data(:)) max(theCircularPSFData.data(:))]);

    visualizedProfile = 'LSF'; % choose between 'midRow' and 'LSF' (integral along Y)
    
    maxRF  = max([...
                max(abs(RF2DData.visualRF(:))) ...
                max(abs(RF2DData.retinalRF(:))) ...
                max(abs(RF2DData.retinalConePoolingRF(:))) ...
                max(abs(RF2DData.visualRFcorrespondingToRetinalConePoolingRF(:)))]);

    switch (visualizedProfile)
        case 'midRow'
            maxProfile = maxRF;

            maxRFcenterConeMap = max([ ...
                max(RF2DData.visualRFcenterConeMap(:)) ...
                max(RF2DData.retinalRFcenterConeMap(:))]);

        case 'LSF'
            maxProfile  = max([...
                max(sum(RF2DData.visualRF,1))...
                max(sum(RF2DData.retinalRF,1)) ...
                max(sum(RF2DData.retinalConePoolingRF,1))...
                max(sum(RF2DData.visualRFcorrespondingToRetinalConePoolingRF,1))]);
        
            maxRFcenterConeMap = max([ ...
                max(sum(RF2DData.visualRFcenterConeMap,1)) ...
                max(sum(RF2DData.retinalRFcenterConeMap,1))]);

        otherwise
            error('Unknown visualized profile: ''%s''.', visualizedProfile);
    end

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [15 15 1350 1000], 'Color', [1 1 1]);

    % The original PSF
    ax = subplot(3,4,1);
    plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID);


    % The circular PSF
    ax = subplot(3,4,2);
    plotPSF(ax, theCircularPSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID);
    title('circularly-averaged PSF')

    % The retinal RF center cone map
    ax = subplot(3,4,3);
    plotRF(ax, thePSFData, RF2DData.retinalRFcenterConeMap, maxRFcenterConeMap, maxRFcenterConeMap, visualizedProfile, maxSpatialSupportDegs, 'retinal RF center');


    % The visual RF center cone map
    ax = subplot(3,4,4);
    plotRF(ax, thePSFData, RF2DData.visualRFcenterConeMap, maxRFcenterConeMap, maxRFcenterConeMap, visualizedProfile, maxSpatialSupportDegs, 'visual RF center');

    
    % The visual RF
    ax = subplot(3,4,5);
    plotRF(ax, thePSFData, RF2DData.visualRF, maxRF, maxProfile, visualizedProfile, maxSpatialSupportDegs, 'visual RF (target)');


    % The retinal RF (via deconvolution of the visual RF
    ax = subplot(3,4,6);
    plotRF(ax, thePSFData, RF2DData.retinalRF, maxRF, maxProfile, visualizedProfile, maxSpatialSupportDegs, 'retinal RF (continuous)');

    % The cone pooling RF 
    ax = subplot(3,4,7);
    plotRF(ax, thePSFData, RF2DData.retinalConePoolingRF, maxRF, maxProfile, visualizedProfile, maxSpatialSupportDegs, 'retinal RF (cone pooling)');


    % The visual RF from the cone pooling retinal RF
    ax = subplot(3,4,8);
    plotRF(ax, thePSFData, RF2DData.visualRFcorrespondingToRetinalConePoolingRF, maxRF, maxProfile, visualizedProfile, maxSpatialSupportDegs, 'visual RF (achieved)');

    % Comparison of target visual and obtained visual RF
    ax = subplot(3,4,9);
    plotProfiles(ax, thePSFData, RF2DData.visualRF, RF2DData.visualRFcorrespondingToRetinalConePoolingRF, visualizedProfile, maxSpatialSupportDegs, 'target vs. achieved visual RF');

   
    % The DoG params
    ax = subplot(3,4,10);
    XYLims = [0.06 6];
    h1 = plot(ax,visualRFparams.RcDegs*60, retinalRFparams.RcDegs*60, ...
        'ro', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], ...
        'LineWidth', 1.0);
    hold(ax, 'on');
    h2 = plot(ax,achievedVisualRFparams.RcDegs*60, retinalRFparams.RcDegs*60, ...
        'ks', 'MarkerSize', 12, 'MarkerFaceColor', [.5 1 .5], ...
        'LineWidth', 1.0);
    
    plot(ax, XYLims, XYLims, 'k-');
    axis(ax, 'square');
    set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
            'XTick', [0.06 0.1 0.3 0.6 1 3 6], 'YTick', [0.06 0.1 0.3 0.6 1 3 6], ...
            'YScale', 'log', 'XScale', 'log', 'FontSize', 14);
    grid(ax, 'on')
    legend(ax,[h1 h2], {'target', 'achieved'}, 'Location', 'NorthOutside', 'NumColumns', 2);
    xlabel(ax, 'visual Rc (arc min)');
    ylabel(ax, 'retinal Rc (arc  min)');

    
    if isfield(retinalRFparams, 'surroundToCenterRcRatio')
        ax = subplot(3,4,11);
        XYLims = [0.5 15];
        h1 = plot(ax,visualRFparams.surroundToCenterRcRatio, ...
            retinalRFparams.surroundToCenterRcRatio, ...
            'ro', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], ...
            'LineWidth', 1.0);
        hold(ax, 'on');
        h2 = plot(ax,achievedVisualRFparams.surroundToCenterRcRatio, ...
            retinalRFparams.surroundToCenterRcRatio, ...
            'ks', 'MarkerSize', 12, 'MarkerFaceColor', [.5 1 .5], ...
            'LineWidth', 1.0);
        plot(ax, XYLims, XYLims, 'k-');
        legend(ax,[h1 h2], {'target', 'achieved'}, 'Location', 'NorthOutside', 'NumColumns', 2);
        
        axis(ax, 'square');
        set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
                'XTick', 0:2:20, 'YTick', 0:2:20, 'FontSize', 14);
        grid(ax, 'on')
        xlabel(ax, 'visual Rs/Rc');
        ylabel(ax, 'retinal Rs/Rc');
    end

    
    if (isfield(retinalRFparams, 'surroundToCenterIntegratedRatio'))
        ax = subplot(3,4,12);
        XYLims = [0.1 10];
    
        h1 = plot(ax,visualRFparams.surroundToCenterIntegratedRatio, ...
            retinalRFparams.surroundToCenterIntegratedRatio, ...
            'ro', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], ...
            'LineWidth', 1.0);
    
        hold(ax, 'on');
        h2 = plot(ax,achievedVisualRFparams.surroundToCenterIntegratedRatio, ...
            retinalRFparams.surroundToCenterIntegratedRatio, ...
            'ks', 'MarkerSize', 12, 'MarkerFaceColor', [.5 1 .5], ...
            'LineWidth', 1.0);
        plot(ax, XYLims, XYLims, 'k-');
        legend(ax,[h1 h2], {'target', 'achieved'}, 'Location', 'NorthOutside', 'NumColumns', 2);
        axis(ax, 'square');
        set(ax, 'XLim', XYLims, 'YLim', XYLims, ...
                'XTick', [0.1 0.3 1 3 10], 'YTick', [0.1 0.3 1 3 10], ...
                'YScale', 'log', 'XScale', 'log', 'FontSize', 14);
        grid(ax, 'on')
        xlabel(ax, 'visual int. S/C ratio');
        ylabel(ax, 'retinal int. S/C ratio');
    end
    

    NicePlot.exportFigToPDF('deconvAnalysis.pdf', hFig, 300);
end

function plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, eccDegs, testSubjectID)
    psfZLevels = 0.05:0.1:0.95;
    contourf(ax,thePSFData.supportX/60, thePSFData.supportY/60, thePSFData.data/maxPSF, psfZLevels);
    hold on;
    midRow = (size(thePSFData.data,1)-1)/2+1;
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.75 + 1.7*thePSFData.data(midRow,:)/maxPSF*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    axis(ax,'image'); axis 'xy';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -0.5:0.05:0.5, 'YTick', -0.5:0.05:0.5, 'CLim', [0 1], 'FontSize', 14);
    grid(ax, 'on');
    xlabel(ax,'degrees');
    title(ax, sprintf('PSF @(%2.0f,%2.0f) degs, subj: %d', eccDegs(1), eccDegs(2), testSubjectID));
    colormap(ax,brewermap(1024, 'greys'));
end


function plotRF(ax, thePSFData, RF, maxRF, maxProfile, visualizedProfile, maxSpatialSupportDegs, titleString)
    %rfZLevels = -0.9:0.1:0.9;
    %contourf(ax,thePSFData.supportX/60, thePSFData.supportY/60, RF/maxRF, rfZLevels);
    imagesc(ax, thePSFData.supportX/60, thePSFData.supportY/60, RF/maxRF);
    hold on;
    switch visualizedProfile
        case 'midRow'
            midRow = (size(RF,1)-1)/2+1;
            theProfile = RF(midRow,:)/maxProfile;
        case 'LSF'
            theProfile = sum(RF,1)/maxProfile;
    end

    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.5 + 1.5*theProfile*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.5 + thePSFData.supportX*0, 'k-');
    axis(ax,'image'); axis 'xy';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -0.5:0.05:0.5, 'YTick', -0.5:0.05:0.5, 'CLim', 0.01*[-1 1], 'FontSize', 14);
    grid(ax, 'on');
    colormap(ax,brewermap(1024, '*RdBu'));
    xlabel(ax,'degrees');
    title(ax, titleString);
end

function plotProfiles(ax, thePSFData, RF1, RF2, visualizedProfile, maxSpatialSupportDegs, titleString)
    
    switch visualizedProfile
        case 'midRow'
            midRow = (size(RF1,1)-1)/2+1;
            theProfile1 = RF1(midRow,:);
            theProfile2 = RF2(midRow,:);
        case 'LSF'
            theProfile1 = sum(RF1,1);
            theProfile2 = sum(RF2,1);
    end
    
    maxProfile = max([max(abs(theProfile1)) max(abs(theProfile2))]);
    theProfile1 = theProfile1 / maxProfile;
    theProfile2 = theProfile2 / maxProfile;

    plot(ax, thePSFData.supportX/60, theProfile1, 'k-', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax, thePSFData.supportX/60, theProfile2, 'r--', 'LineWidth', 1.5);
    plot(ax, thePSFData.supportX/60, theProfile1-theProfile2, 'b-', 'LineWidth', 1.0);
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', [-0.5 1], ...
            'XTick', -0.5:0.05:0.5, 'YTick', -0.6:0.1:1, 'FontSize', 14);
    axis(ax, 'square');
    grid(ax, 'on');
    legend({'target', 'achieved', 'residual'});
    minV = max(abs(RF1(:)))*0.01;
    idx = find(abs(RF1(:))>minV);
    RMSE = 100*sqrt(mean(((RF1(idx)-RF2(idx))./RF1(idx)).^2));
    title(sprintf('RMSE: %2.1f%%',RMSE));
    xlabel(ax,'degrees');
end