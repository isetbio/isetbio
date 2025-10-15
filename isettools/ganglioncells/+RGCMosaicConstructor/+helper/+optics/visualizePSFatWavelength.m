function visualizePSFatWavelength(ax, thePSF, theVisualizedWavelengthIndex, maxPSF, theTitle, varargin)

	% Parse input
    p = inputParser;
    % Optional params
    p.addParameter('XLimsArcMin', [-10 10], @isnumeric);
    p.addParameter('YLimsArcMin', [-10 10], @isnumeric);
    p.addParameter('XTicksArcMin', -20:5:20, @isnumeric);
    p.addParameter('YTicksArcMin', -20:5:20, @isnumeric);
    p.parse(varargin{:});
    XLimsArcMin = p.Results.XLimsArcMin;
    YLimsArcMin = p.Results.YLimsArcMin;
    XTicksArcMin = p.Results.XTicksArcMin;
    YTicksArcMin = p.Results.YTicksArcMin;

	theVisualizedPSFslice = squeeze(thePSF.data(:,:,theVisualizedWavelengthIndex));
	imagesc(ax, thePSF.supportX, thePSF.supportY, theVisualizedPSFslice / maxPSF);
	hold(ax, 'on');
	plot(ax, XLimsArcMin, [0 0], 'r-');
	plot(ax, [0 0], YLimsArcMin, 'r-');
	axis(ax, 'square'); axis(ax, 'image');
	set(ax, 'CLim', [0 1], 'XLim', XLimsArcMin, 'YLim', YLimsArcMin, 'XTick', XTicksArcMin, 'YTick', YTicksArcMin, 'FontSize', 16);
	xlabel(ax, 'space, x (arc min)');
	ylabel(ax, 'space, y (arc min)');
	grid(ax, 'on');
    colormap(ax, 'gray');

	if (~isempty(theTitle))
	    title(ax, theTitle);
	else
		title(sprintf('max of PSF at visualized wavelength: %f', max(thePSFsliceAtWavelengthIndexOfPeakAmplitude(:))));
	end
	drawnow;
end

