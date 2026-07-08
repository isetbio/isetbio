function visualizeSceneRadiance(spatialSupport, spatialSupportUnits, ...
    photons, wavelengthSupport, wavelengthBandsToVisualize, varargin)

p = inputParser;
p.addParameter('contrastProfilesOnly', false, @islogical);
% Parse input
p.parse(varargin{:});
contrastProfilesOnly = p.Results.contrastProfilesOnly;

% Method to visualize the scene radiance at certain wavelength bands
spatialSupportX = squeeze(spatialSupport(1,:,1));
spatialSupportY = squeeze(spatialSupport(:,1,2));
photonRange = [min(photons(:)) max(photons(:))];

subplotRows = 3;
subplotCols = 3;
subplotPos = NicePlot.getSubPlotPosVectors(...
    'rowsNum', subplotRows, 'colsNum', subplotCols, ...
    'heightMargin',  0.07, 'widthMargin',    0.01, ...
    'leftMargin',     0.01, 'rightMargin',    0.01, ...
    'bottomMargin',   0.05, 'topMargin',      0.04);

figure(); clf;
cmap = brewermap(1024,'*Spectral');
colormap(cmap);
hProgressBar = waitbar(0,'Getting photon rates ...');
for iBand = 1:numel(wavelengthBandsToVisualize)
    if (iBand > subplotRows*subplotCols)
        continue
    end
    waitbar(iBand/numel(wavelengthBandsToVisualize), hProgressBar, ...
        sprintf('Getting photon rates at %2.0f nm', wavelengthBandsToVisualize(iBand)));
    row = floor((iBand-1)/subplotCols)+1; col = mod(iBand-1,subplotCols)+1;
    subplot('Position', subplotPos(row,col).v);
    % find index of spectral slices requested
    [~,visualizedWavelengthIndex] = ...
        min(abs(wavelengthSupport-wavelengthBandsToVisualize(iBand)));
    
    if (contrastProfilesOnly)
        midRow = floor(size(photons,1)/2)+1;
        profile = squeeze(photons(midRow,:,visualizedWavelengthIndex));
        % backgroundPhotonstaken from the origin (1,1)
        backgroundPhotons = photons(1,1,visualizedWavelengthIndex);
        WeberContrastProfile = (profile-backgroundPhotons)/backgroundPhotons;
        plot(spatialSupportX, WeberContrastProfile, 'r-', 'LineWidth', 1.5);
        grid on
        set(gca, 'XTick', [spatialSupportX(1) 0 spatialSupportX(end)], 'YTick', [-1:0.5:1]);
        set(gca, 'YLim', [-1 1], 'XLim', [spatialSupportX(1) spatialSupportX(end)]);
        xtickformat('%0.2f');
        if (row < subplotRows)
            set(gca, 'XTickLabel', {});
        end
        axis 'square'
    else
        imagesc(spatialSupportX, spatialSupportY, ...
            squeeze(photons(:,:,visualizedWavelengthIndex)), photonRange);
        axis 'image';
        set(gca, 'XTick', [], 'YTick', []);
        colorbar();
    end

    set(gca, 'FontSize', 16);
    title(sprintf('%2.0fnm', wavelengthSupport(visualizedWavelengthIndex))); 
end
close(hProgressBar);
drawnow;

end