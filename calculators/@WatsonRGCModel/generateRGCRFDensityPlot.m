function generateRGCRFDensityPlot(obj, RGCRFDensityFunctionHandle, eccDegs)
    % Plot the total RGC RF density along the temporal meridian
    plot(eccDegs, RGCRFDensityFunctionHandle(eccDegs, 'temporal meridian', 'RFs per deg2'), ...
        'r-', 'LineWidth', obj.figurePrefs.lineWidth); hold on;
    plot(eccDegs, RGCRFDensityFunctionHandle(eccDegs, 'superior meridian', 'RFs per deg2'), ...
        'b-', 'LineWidth', obj.figurePrefs.lineWidth);
    plot(eccDegs, RGCRFDensityFunctionHandle(eccDegs, 'nasal meridian', 'RFs per deg2'), ...
        'g-', 'LineWidth', obj.figurePrefs.lineWidth);
    plot(eccDegs, RGCRFDensityFunctionHandle(eccDegs, 'inferior meridian', 'RFs per deg2'), ...
        'k-', 'LineWidth', obj.figurePrefs.lineWidth);
    legend({'temporal', 'superior', 'nasal', 'inferior'}, 'Location', 'SouthWest');
    xlabel('eccentricity (degs)', 'FontAngle', obj.figurePrefs.fontAngle);
    
    if (isequal(func2str(RGCRFDensityFunctionHandle),'@(varargin)obj.midgetRGCRFDensity(varargin{:})'))
        ylabel('midget RGC density (# of RFs/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
        title('Midget RGC RF density (ON+OFF) as a function of eccentricity for four meridians');
    elseif (isequal(func2str(RGCRFDensityFunctionHandle),'@(varargin)obj.totalRGCRFDensity(varargin{:})'))
        ylabel('total RGC density (# of RFs/deg^2)', 'FontAngle', obj.figurePrefs.fontAngle);
        title('Total RGC RF density as a function of eccentricity for four meridians');
    else
        error('RGCRFDensityFunctionHandle must be either ''@obj.midgetRGCRFDensity'' or ''@obj.totalRGCRFDensity''.');
    end
    set(gca, 'XLim', [0.07 120], 'YLim', [1 35*1000], ...
        'XTick', [0.1 0.5 1 5 10 50 100], 'YTick', [1 10 100 1000 10000], ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', obj.figurePrefs.fontSize);
    grid(gca, obj.figurePrefs.grid);
end


