function [udata, g] = oiPlot(oi, pType, roiLocs, varargin)
% Gateway routine for plotting optical image (oi) properties
%
% Syntax:
%   [udata, g] = oiPlot([oi], [pType], [roiLocs], [varargin])
%
% Description:
%    Gateway routine to plot the irradiance or illuminance data in the
%    optical image. There are many options.
%
%    The data shown in the plot are generally returned in udata. The data
%    can also be retrieved from the figure itself, using the call
%    uData = get(figHandle, 'userdata');
%
%    Inputs are the optical image (oi), the plot type (pType), in some
%    cases a position or ROI locations is required (xy) and in some cases
%    other arguments can be included to make plotting from scripts possible
%    without requiring user intervention (e.g., grid spacing in irradiance
%    image with grid).
%
% Inputs:
%    oi       - (Optional) The optical Image structure. Default is to
%               retrieve an oi object using vcGetObject.
%    pType    - (Optional). The plot type. Default is 'illuminance hline'.
%               There are a vast number of types available, listed below,
%               sorted by category:
%           Irradiance (Irr)
%               {'irradiance photons roi'} - Irr within an ROI of the image
%               {'irradiance energy roi'} - Irr within an ROI of the image
%               {'irradiance vline'} - Vertical line spectral irradiance
%                                      (photons) - (space x wavelength)
%               {'irradiance hline'} - Horizontal line spectral irradiance
%                                      (photons) - (space x wavelength)
%               {'irradiance fft'} - 2D FFT of radiance at some wavelength
%               {'irradiance image grid'} - Show spatial grid on irr image
%               {'irradiance image no grid'} - Show irr image without grid
%               {'irradiance waveband image'} - Irr image within a band
%           Illuminance (Ill)
%               {'illuminance mesh log'} - Mesh plot of image log ill
%               {'illuminance mesh linear'} - Mesh plot of image ill
%               {'illuminance fft'} - 2D FFT of ill
%               {'illuminance hline'} - Horizontal line ill
%               {'illuminance fft hline'} - Horizontal line ill fft
%               {'illuminance vline'} - Vertical line luminance
%               {'illuminance fft vline'} - Vertical line luminance FFT
%               {'illuminance roi'} - Histogram of ill in an ROI
%           CIE
%               {'chromaticity roi'} - CIE xy in a region of interest
%           Contrast (Con)
%               {'contrast hline'} - Horizontal line con at a wavelength
%               {'contrast vline'} - Vertical line con
%           Depth
%               {'depth map'} - If has a depth map, plot as a mesh
%               {'depth map contour'} - If has a depth map, plot as a mesh
%           Optics related:
%               {'otf'} - Optical transfer function, units are lines/mm
%               {'otf 550'} - OTF at 550 nm
%               {'otf wavelength'} - One dimensional cut through the OTF at
%                                    all wavelengths. Units are cycles/mm
%               {'psf'} - Point spread function at selected wavelength
%               {'psf 550'} - PSF at 550nm spatial units are microns
%               {'ls wavelength'} - Line spread function at all wavelengths
%                   Peak spatial frequency can be set for the OTF (default:
%                   3 * incoherent cutoff). Number of spatial samples to
%                   plot in the line spread can be set (default: 40).
%               {'lens transmittance'} - Spectral lens transmittance.
%                   Computed from the lens density in the human case.
%    roiLocs  - (Optional) Region of Interest Locations. Default depends on
%               the plot type in order select the region of interest, which
%               is sometimes a line and sometimes a rectangle.
%    varargin - Additional arguments passed in, such as:
%           {wave, gSpacing} - Specify wavelength and/or grid spacing.
%
% Outputs:
%    udata   - User data structure
%    g       - Figure/Graph Handle
%
% Optional key/value pairs:
%    Needs to be populated.
%
% Notes:
%    * [Note: JNM - roiLocs are listed as optional, however a number of the
%      functions will break if they have not been provided because there is
%      not language in place to generate them for all of the cases.  What
%      needs to be done in such cases is to add the pType to one of the two
%      cases at the top of the routine, where the users is asked to select
%      and roi.  One case is when the roi should be a line, the other is
%      when it should be a rectangle.
%    * [Note: DHB - Someday might convert code that gets and plots lsf to
%       use PTB wrapper routines in external, to improve overall
%       consistency.
%
% See Also:
%    t_oiPlot, scenePlot
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/11/17  jnm  Formatting
%    12/22/17  dhb  Use opticsGet(...'diffractionlimitedpsfdata'...) to get
%                   diffraction limited psf, not opticsGet(...'psf'...).
%                   The latter seemed unfortunately named. (Reverted,
%                   BW).
%    12/28/17  dhb  Separated out into separate grouped switch statements.
%    12/30/17  dhb  Went through and verified that various OTF things are
%                   done in a manner consistent with recent changes to
%                   optics, and tried to comment key points more fully.
%                   Removed note that I should take a look.
%              dhb  Add ifftshift for fft calls, since we view data as an
%                   image with zero in the center.  Doesn't matter because
%                   only the amplitude is being plotted, but seemed like
%                   good coding practice.
%    01/24/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip.  The tutorial is auto run when we test tutorials.
    t_oiPlot
%}

if notDefined('oi'), oi = vcGetObject('OI'); end
if notDefined('pType'), pType = 'hlineilluminance'; end

% Reformat the parameter - no spaces, all lower case
pType = ieParamFormat(pType);

% Set spatial locations to plot, if not set already
if notDefined('roiLocs')
    % oiWindow;
    switch pType
        case {  'irradiancevline', 'vline', 'vlineirradiance', ...
                'irradiancehline', 'hline', 'hlineirradiance' , ...
                'illuminancehline', 'horizontallineilluminance', ...
                'hlineilluminance', ...
                'illuminanceffthline', 'illuminancefftvline', ...
                'illuminancevline', 'vlineilluminance', ...
                'contrasthline', 'hlinecontrast', ...
                'contrastvline', 'vlinecontrast'}
            roiLocs = vcLineSelect(oi);

        case {'irradianceenergyroi', 'irradiancephotonsroi', ...
                'chromaticityroi', 'illuminanceroi'}
            roiLocs = vcROISelect(oi);

        otherwise
            % There are cases that don't need a position
    end
end

% Make the plot window and use this default gray scale map.
g = vcNewGraphWin;
mp = 0.4 * gray + 0.4 * ones(size(gray));
colormap(mp);

% Here we have a series of switches, each grouping together
% plots of different basic types
%
% The first group is irradiance related
isIrradiancePlot = true;
switch pType
    case {'irradiancephotonsroi'}
        %[uData, g] = oiPlot(oi, 'irradiance photons roi', roiLocs);
        udata = oiPlotIrradiance(oi, 'photons', roiLocs);

    case {'irradianceenergyroi'}
        %[uData, g] = oiPlot(oi, 'irradiance energy roi', roiLocs);
        udata = oiPlotIrradiance(oi, 'energy', roiLocs);

    case {'irradiancehline', 'hline', 'hlineirradiance'}
        % oiPlot('irradiance hline')
        data = oiGet(oi, 'photons');
        if isempty(data)
            warndlg(sprintf('Photon data are unavailable.'));
            return;
        end

        wave = oiGet(oi, 'wave');
        data = squeeze(data(roiLocs(2), :, :));
        if isa(data, 'single'), data = double(data); end

        posMicrons = oiSpatialSupport(oi, 'um');

        if size(data, 1) == 1
            % Manage monochrome data
            plot(posMicrons.x, data');
            xlabel('Position (um)');
            ylabel('Irradiance (q/s/nm/m^2)');
            title('Monochrome image')
            grid on;
            set(gca, 'xtick', ieChooseTickMarks(posMicrons.x))
        else
            mesh(posMicrons.x, wave, double(double(data')));
            xlabel('Position (um)');
            ylabel('Wavelength (nm)');
            zlabel('Irradiance (q/s/nm/m^2)');
            grid on;
            set(gca, 'xtick', ieChooseTickMarks(posMicrons.x))
        end

        udata.wave = wave;
        udata.pos = posMicrons.x;
        udata.data = double(data');
        udata.cmd = 'mesh(pos, wave, data)';
        set(g, 'Name', sprintf('ISET GraphWin: line %.0f', roiLocs(2)));
        colormap(jet)

    case {'irradiancevline', 'vline', 'vlineirradiance', }
        % oiPlot(oi, 'irradiance vline')
        data = oiGet(oi, 'photons');
        if isempty(data)
            warndlg(sprintf('Photon data are unavailable.'));
            return;
        end

        wave = oiGet(oi, 'wave');
        data = squeeze(data(:, roiLocs(1), :));
        if isa(data, 'single'), data = double(data); end

        posMicrons = oiSpatialSupport(oi, 'microns');

        if size(data, 2) == 1
            plot(posMicrons.y, data);
            xlabel('Position (mm)');
            ylabel('Irradiance (q/s/nm/m^2)');
            title('Monochrome image')
            grid on;
            set(gca, 'xtick', ieChooseTickMarks(posMicrons.y))
        else
            mesh(posMicrons.y, wave, double(data'));
            xlabel('Position (um)');
            ylabel('Wavelength (nm)');
            zlabel('Irradiance (q/s/nm/m^2)');
            grid on;
            set(gca, 'xtick', ieChooseTickMarks(posMicrons.y))
        end

        % Attach data to the figure
        udata.wave = wave;
        udata.pos = posMicrons.y;
        udata.data = double(data');
        set(g, 'Name', sprintf('Line %.0f', roiLocs(1)));

        colormap(jet)

    case {'irradiancefft'}
        % plot(oi, 'irradiance fft', roiLocs, wave)
        % This is the fft of the region at the selected wavelength.
        %
        % Default: roiLocs - whole image, roiLocs not tested adequately.
        %          wave    - middle wavelength
        %
        % The mean is not included in the graph to help with the dynamic
        % range.
        %
        % Axis range could be better.
        if isempty(varargin)
            wave = oiGet(oi, 'wave');
            selectedWave = wave(round(length(wave) / 2));
        else
            selectedWave = varargin{1};
        end

        data = oiGet(oi, 'photons', selectedWave);
        if isempty(data)
            warndlg(sprintf('Photon data are unavailable.'));
            return;
        end

        if isa(data, 'single'), data = double(data); end
        sz = size(data);

        % Remove the mean
        data = data - mean(data(:));

        % Plot and attach data to figure.
        %
        % The ifftshift in front of the fft2 call is because we view what
        % is in data as an image with (0,0) at the center.  It doesn't
        % matter since we're not looking at the phase of the fft, but seems
        % best to have as good coding practice.
        udata.x = 1:sz(2);
        udata.y = 1:sz(1);
        udata.z = fftshift(abs(fft2(ifftshift(data))));
        udata.cmd = 'mesh(x, y, z)';
        mesh(udata.x, udata.y, udata.z);
        xlabel('Cycles/ROI-image');
        ylabel('Cycles/ROI-image');
        zlabel('Amplitude');
        str = sprintf('Amplitude spectrum at %.0f nm', selectedWave);
        title(str);
        set(g, 'Name', sprintf('Irradiance with grid'));
        colormap(jet)

    case {'irradiancewavebandimage', 'irradiancewavebandimagegrid', ...
            'irradianceimagewave', 'irradianceimagewavegrid'}
        % oiPlot(oi, 'irradianceImageWave', wave, gSpacing);
        if isempty(varargin)
            wave = 500;
        else
            wave = varargin{1};
        end

        irrad = oiGet(oi, 'photons', wave);
        sz = oiGet(oi, 'size');
        spacing = oiGet(oi, 'sampleSpacing', 'um');

        % This is probably now a spatial support oiGet ...
        xCoords = spacing(2) * (1:sz(2));
        xCoords = xCoords - mean(xCoords);
        yCoords = spacing(1) * (1:sz(1));
        yCoords = yCoords - mean(yCoords);
        suggestedSpacing = round(max(xCoords(:)) / 5);
        if length(varargin) == 2
            gSpacing = varargin{2};
        else
            gSpacing = ieReadNumber('Enter grid spacing (um)', ...
                suggestedSpacing, '%.2f');
            if isempty(gSpacing), return; end
        end

        rgb = imageSPD(irrad, wave);
        imagesc(xCoords, yCoords, rgb);
        xlabel('Position (um)');
        ylabel('Position (um)');

        udata.irrad = irrad;
        udata.xCoords = xCoords;
        udata.yCoords = yCoords;

        xGrid = (0:gSpacing:round(max(xCoords)));
        tmp = -1 * fliplr(xGrid);
        xGrid = [tmp(1:(end - 1)), xGrid];
        yGrid = (0:gSpacing:round(max(yCoords)));
        tmp = -1 * fliplr(yGrid);
        yGrid = [tmp(1:(end - 1)), yGrid];

        set(gca, 'xcolor', [.5 .5 .5]);
        set(gca, 'ycolor', [.5 .5 .5]);
        set(gca, 'xtick', xGrid, 'ytick', yGrid);
        grid on
        set(g, 'Name', sprintf('Image with grid'));

    case {'irradianceimagegrid', 'irradianceimagewithgrid', ...
            'irradianceimage'}
        % oiPlot(oi, 'irradianceImage', sampleSpacing-um);
        irrad = oiGet(oi, 'photons');
        wave = oiGet(oi, 'wave');
        sz = oiGet(oi, 'size');
        spacing = oiGet(oi, 'sampleSpacing', 'um');

        % This is probably now a spatial support oiGet ...
        xCoords = spacing(2) * (1:sz(2));
        xCoords = xCoords - mean(xCoords);
        yCoords = spacing(1) * (1:sz(1));
        yCoords = yCoords - mean(yCoords);
        if length(varargin) >= 1
            gSpacing = varargin{1};
        else
            suggestedSpacing = round(max(xCoords(:)) / 5);
            gSpacing = ieReadNumber('Enter grid spacing (um)', ...
                suggestedSpacing, '%.2f');
            if isempty(gSpacing), return; end
        end

        nWave = oiGet(oi, 'nwave');
        wList = oiGet(oi, 'wavelength');
        [row, col] = size(irrad);
        if nWave > 1
            imageSPD(irrad, wave, [], row, col, 1, xCoords, yCoords);
        else
            imageSPD(irrad, wList, [], row, col, 1, xCoords, yCoords);
        end
        xlabel('Position (um)');
        ylabel('Position (um)');

        udata.irrad = irrad;
        udata.xCoords = xCoords;
        udata.yCoords = yCoords;

        xGrid = (0:gSpacing:round(max(xCoords)));
        tmp = -1 * fliplr(xGrid);
        xGrid = [tmp(1:(end - 1)), xGrid];
        yGrid = (0:gSpacing:round(max(yCoords)));
        tmp = -1 * fliplr(yGrid);
        yGrid = [tmp(1:(end - 1)), yGrid];

        set(gca, 'xcolor', [.5 .5 .5]);
        set(gca, 'ycolor', [.5 .5 .5]);
        set(gca, 'xtick', xGrid, 'ytick', yGrid);
        grid on
        set(g, 'Name', sprintf('Irradiance with grid'));

    case {'irradianceimagenogrid'}
        % oiPlot(oi, 'irradianceImage', sampleSpacing-um);
        irrad = oiGet(oi, 'photons');
        wave = oiGet(oi, 'wave');
        sz = oiGet(oi, 'size');
        spacing = oiGet(oi, 'sampleSpacing', 'um');

        % This is probably now a spatial support oiGet ...
        xCoords = spacing(2) * (1:sz(2));
        xCoords = xCoords - mean(xCoords);
        yCoords = spacing(1) * (1:sz(1));
        yCoords = yCoords - mean(yCoords);

        nWave = oiGet(oi, 'nwave');
        wList = oiGet(oi, 'wavelength');
        [row, col] = size(irrad);
        if nWave > 1
            imageSPD(irrad, wave, [], row, col, 1, xCoords, yCoords);
        else
            imageSPD(irrad, wList, [], row, col, 1, xCoords, yCoords);
        end
        xlabel('Position (um)');
        ylabel('Position (um)');
        grid off

        udata.irrad = irrad;
        udata.xCoords = xCoords;
        udata.yCoords = yCoords;

    otherwise
        % If we're here, pType did not specify an irradiance plot
        isIrradiancePlot = false;
end

% Switch for lluminance and chromaticity plots
isIlluminancePlot = true;
switch (pType)
    case {'illuminanceroi'}
        % Histogram of illuminance in an ROI
        udata = oiPlotCIE(oi, 'illuminance', roiLocs);

    case {'chromaticityroi'}
        % Graph of chromaticity coords in an ROI
        udata = oiPlotCIE(oi, 'chromaticity', roiLocs);

    case {'illuminancehline', 'horizontallineilluminance', ...
            'hlineilluminance'}
        % oiPlot(oi, 'illuminance hline')
        data = oiGet(oi, 'illuminance');
        if isempty(data)
            warndlg(sprintf('Illuminance data are unavailable.'));
            return;
        end
        illum = data(roiLocs(2), :);
        posMicrons = oiSpatialSupport(oi, 'um');

        plot(posMicrons.x, illum);
        xlabel('Position (um)');
        ylabel('Illuminance (lux)');
        grid on;
        set(gca, 'xtick', ieChooseTickMarks(posMicrons.x))

        udata.pos = posMicrons.x;
        udata.data = illum';
        udata.cmd = 'plot(pos, illum)';
        set(g, 'Name', sprintf('Line %.0f', roiLocs(2)));

    case {'illuminancemeshlog'}
        % Mesh plot of image log illuminance
        udata = oiPlotIlluminanceMesh(oi, 'log');

    case {'illuminancemeshlinear'}
        % Mesh plot of image illuminance
        udata = oiPlotIlluminanceMesh(oi, 'linear');

    case {'illuminanceffthline'}
        % oiPlot(oi, 'illuminance fft hline')
        % The mean is removed to keep the dynamic range reasonable.

        data = oiGet(oi, 'illuminance');
        if isempty(data)
            warndlg(sprintf('Illuminance data are unavailable.'));
            return;
        end
        illum = data(roiLocs(2), :);
        illum = illum - mean(illum(:));
        xPosMM = oiSpatialSupport(oi, 'mm');

        % Compute amplitude spectrum in units of millimeters
        normalize = 1;
        [freq, fftIllum] = ieSpace2Amp(xPosMM.x, illum, normalize);

        plot(freq, fftIllum, 'r-');
        xlabel('Cycles/mm');
        ylabel('Normalized amplitude');
        grid on

        udata.freq = freq;
        udata.data = fftIllum;
        udata.cmd = 'plot(freq, data, ''r-'')';
        set(g, 'Name', sprintf('Line %.0f', roiLocs(2)));

    case {'illuminancevline', 'vlineilluminance'}
        % oiPlot(oi, 'illuminance vline')

        data = oiGet(oi, 'illuminance');
        if isempty(data)
            warndlg(sprintf('Illuminance data are unavailable.'));
            return;
        end
        illum = data(:, roiLocs(1));
        posMicrons = oiSpatialSupport(oi, 'um');

        plot(posMicrons.y, illum);
        xlabel('Position (um)');
        ylabel('Illuminance (lux)');
        grid on;
        set(gca, 'xtick', ieChooseTickMarks(posMicrons.y))

        udata.pos = posMicrons.y;
        udata.data = illum';
        udata.cmd = 'plot(pos, illum)';
        set(g, 'Name', sprintf('Line %.0f', roiLocs(1)));

    case {'illuminancefftvline'}
        % oiPlot(oi, 'illuminance fft vline')
        % space = oiGet(oi, 'spatialSupport');

        data = oiGet(oi, 'illuminance');
        if isempty(data)
            warndlg(sprintf('Illuminance data are unavailable.'));
            return;
        end
        illum = data(:, roiLocs(1));
        yPosMM = oiSpatialSupport(oi, 'mm');

        % Compute amplitude spectrum in units of millimeters
        normalize = 1;
        [freq, fftIllum] = ieSpace2Amp(yPosMM.y, illum, normalize);

        plot(freq, fftIllum, 'r-');
        xlabel('Cycles/mm');
        ylabel('Normalized amplitude');
        grid on

        udata.freq = freq;
        udata.data = fftIllum;
        udata.cmd = 'plot(freq, data, ''r-'')';
        set(gcf, 'Name', sprintf('Line %.0f', roiLocs(1)));

    case {'illuminancefft', 'fftilluminance'}
        % oiPlot(oi, 'illuminance fft')

        % The ifftshift in front of the fft2 call is because we view what
        % is in data as an image with (0,0) at the center.  It doesn't
        % matter since we're not looking at the phase of the fft, but seems
        % best to have as good coding practice.
        data = oiGet(oi, 'illuminance');
        sz = size(data);
        udata.x = 1:sz(2);
        udata.y = 1:sz(1);
        udata.z = fftshift(abs(fft2(ifftshift(data))));
        udata.cmd = 'mesh(x, y, z)';
        mesh(udata.x, udata.y, udata.z);
        xlabel('Cycles/image');
        ylabel('Cycles/image');
        zlabel('Amplitude');
        title('Illuminance amplitude spectrum');

    otherwise
        % If w're here, it wasn't an illuminance related plot
        isIlluminancePlot = false;
end

% Contrast related plots
isContrastPlot = true;
switch (pType)
    case {'contrasthline', 'hlinecontrast'}
        % oiPlot(oi, 'contrast hline')
        % Plot percent contrast (difference from the mean as a percentage
        % of the mean).

        data = oiGet(oi, 'photons');
        if isempty(data)
            warndlg(sprintf('Photon data are unavailable.'));
            return;
        end
        data = squeeze(data(roiLocs(2), :, :));
        if isa(data, 'single'), data = double(data); end

        % Percent contrast
        mn = mean(data(:));
        if mn == 0
            warndlg('Zero mean. Cannot compute contrast.');
            return;
        end
        data = 100 * (data - mn) / mn;

        posMicrons = oiSpatialSupport(oi, 'microns');

        wave = oiGet(oi, 'wave');

        mesh(posMicrons.x, wave, double(data'));
        xlabel('Position (um)');
        ylabel('Wavelength (nm)');
        zlabel('Percent contrast');
        grid on;
        set(gca, 'xtick', ieChooseTickMarks(posMicrons.x))
        udata.wave = wave;
        udata.pos = posMicrons.x;
        udata.data = double(data');
        udata.cmd = 'mesh(pos, wave, data)';
        set(g, 'Name', sprintf('Line %.0f', roiLocs(2)));
        colormap(jet)

    case {'contrastvline', 'vlinecontrast'} % Done
        % oiPlot(oi, 'contrast vline')
        % Plot percent contrast (difference from the mean as a percentage
        % of the mean).
        data = oiGet(oi, 'photons');
        if isempty(data)
            warndlg(sprintf('Photon data are unavailable.'));
            return;
        end

        wave = oiGet(oi, 'wave');
        data = squeeze(data(:, roiLocs(1), :));
        if isa(data, 'single'), data = double(data); end

        % Percent contrast
        mn = mean(data(:));
        if mn == 0
            warndlg('Zero mean. Cannot compute contrast.');
            return;
        end
        data = 100 * (data - mn) / mn;

        posMicrons = oiSpatialSupport(oi, 'microns');

        mesh(posMicrons.y, wave, double(data'));
        xlabel('Position (um)');
        ylabel('Wavelength (nm)');
        zlabel('Irradiance (q/s/nm/m^2)')
        zlabel('Percent contrast')
        grid on;
        set(gca, 'xtick', ieChooseTickMarks(posMicrons.y))

        udata.wave = wave;
        udata.pos = posMicrons.y;
        udata.data = double(data');
        udata.cmd = 'mesh(pos, wave, data)';
        set(g, 'Name', sprintf('Line %.0f', roiLocs(1)));

    otherwise
        isContrastPlot = false;
end

% Depth related
isDepthPlot = true;
switch (pType)
    case {'depthmap'}
        % oiPlot(oi, 'depth map')
        dmap = oiGet(oi, 'depth map');
        if isempty(dmap)
            close(g);
            error('No depth map')
        else
            imagesc(dmap); colormap(flipud(gray))
            namestr = sprintf('Depth map (max = %.1f)', max(dmap(:)));
            axis off; axis image; set(g, 'Name', namestr);
        end
        colorbar;
        udata.dmap = dmap;

    case {'depthmapcontour', 'depthcontour'}
        % oiPlot(oi, 'depth contour')
        dmap = oiGet(oi, 'depth map');
        dmap = ieScale(dmap, 0, 1);
        mx = max(dmap(:));
        drgb = cat(3, dmap, dmap, dmap);

        image(drgb);
        colormap(flipud(gray));
        hold on
        n = 4;
        v = (1:n) / n;
        contour(dmap, v);
        hold off
        namestr = sprintf('ISET: Depth map (max = %.1f m)', mx);
        axis off;
        set(g, 'Name', namestr);

    otherwise
        isDepthPlot = false;
end

% Optics related
isOpticsPlot = true;
switch (pType)
    case {'lenstransmittance'}
        % oiPlot(oi, 'lens transmittance')
        % If human, uses the lens object attached to oi.
        w = oiGet(oi, 'wavelength');
        set(g, 'Name', 'ISETBIO:  Lens');

        if checkfields(oi, 'optics', 'lens')
            % human case
            t = oiGet(oi, 'lens transmittance');

            plot(w, t);
            grid on
            xlabel('Wavelength (nm)');
            ylabel('Transmittance');
            d = oiGet(oi, 'lens density');
            title(sprintf('Lens transmittance (density %.2f)', d))
        else
            %Diffraction case
            if isempty(w)
                % No calculation yet, so put in a dummy wavelength array
                % just for this plot. It is not saved.
                w = 400:10:700;
                oi = oiSet(oi, 'optics wave', w);
            end
            t = oiGet(oi, 'optics transmittance');
            plot(w, t);
            grid on
            xlabel('Wavelength (nm)');
            ylabel('Transmittance');
        end

        udata.wavelength = w;
        udata.transmittance = t;
    case {'otf', 'otfanywave'}
        % User asked to select a wavelength
        % Optical transfer function, units are lines/mm
        % oiPlot(oi, 'otf', [], 420);
        optics = oiGet(oi, 'optics');
        opticsModel = opticsGet(optics, 'model');
        switch lower(opticsModel)
            case 'raytrace'
                rtPlot(oi, 'otf');
            otherwise
                if isempty(varargin)
                    udata = oiPlotOTF(oi, 'otf');
                else
                    w = varargin{1};
                    udata = oiPlotOTF(oi, 'otf', w);
                end
        end
        set(g, 'userdata', udata);
        set(g, 'name', 'OTF');
        colormap(jet)

    case {'otf550'}
        % OTF at 550 nm
        udata = oiPlotOTF(oi, 'otf 550');
        set(g, 'userdata', udata);
        set(g, 'name', 'OTF 550');
        colormap(jet)

    case {'psf'}
        % Point spread function at selected wavelength
        % oiPlot(oi, 'psf', [], 420);
        if isempty(varargin)
            udata = oiPlotOTF(oi, 'psf');
        else
            w = varargin{1};
            udata = oiPlotOTF(oi, 'psf', w);
        end
        set(g, 'userdata', udata);
        namestr = sprintf('ISET: %s', oiGet(oi, 'name'));
        set(g, 'Name', namestr);
        colormap(jet)

    case {'psf550'}
        % PSF at 550nm spatial units are microns
        udata = oiPlotOTF(oi, 'psf 550');
        set(g, 'userdata', udata);
        namestr = sprintf('ISET: %s', oiGet(oi, 'name'));
        set(g, 'Name', namestr);
        colormap(jet)

    case {'lswavelength', 'lsfwavelength'}
        % uData = oiPlot(oi, pType, [], nSpatialSamps)
        % the nSpatialSamps part isn't working.
        %
        % Line spread function at all wavelengths.
        %    Peak spatial frequency can be set for the OTF (default: 3 *
        %    incoherent cutoff). Number of spatial samples to plot in the
        %    line spread can be set (default: 40).
        if ~isempty(varargin), nSamps = varargin{1};
        else,                  nSamps = 40;
        end
        udata = oiPlotOTF(oi, 'ls wavelength', [], nSamps);
        set(g, 'userdata', udata);
        set(g, 'name', 'LS by Wave');
        colormap(jet)

    case{'otfwavelength', 'mtfwavelength'}
        % One dimensional otf at all wavelengths as  mesh plot.
        % Units are cycles/mm
        optics = oiGet(oi, 'optics');
        opticsModel = opticsGet(optics, 'model');
        switch lower(opticsModel)
            case 'raytrace'
                % Not what the user asked for. Must fix. Add varargin,
                % and make the right plot. This isn't it.
                rtPlot(oi, 'otf');
            otherwise
                udata = oiPlotOTF(oi, 'otf wavelength');
                set(g, 'userdata', udata);
        end
        set(g, 'name', 'OTF by Wave');
        colormap(jet)

    otherwise
        isOpticsPlot = false;
end

% Check that a known pType was specified
if (~isIrradiancePlot && ~isIlluminancePlot && ~isContrastPlot && ...
        ~isDepthPlot && ~isOpticsPlot)
            error('Unknown oiPlot type %s.', pType);
end

% Set user data for return
if exist('udata', 'var'), set(gcf, 'userdata', udata); end

end

function udata = oiPlotIrradiance(oi, dataType, roiLocs)
% Plot mean irradiance within a selected ROI of the optical image window
%
% Syntax:
%   udata = oiPlotIrradiance(oi, [dataType], roiLocs)
%
% Description:
%    Plot the average optical image irradiance within a selected ROI. The
%    default data type is photons. If the optical image is not  monochrome,
%    the results are plotted in the GRAPHWIN. Otherwise, the mean
%    irradiance within the ROI is displayed in a message box.
%
% Inputs:
%    oi       - The optical image
%    dataType - (Optional) The data type. Default is photons. The options
%               are 'photons' and 'energy'.
%    roiLocs  - Region of Interest Locations.
%
% Outputs:
%    udata    - User Data structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/12/17  jnm  Formatting
%    12/25/17   BW  Fixed by adding roiLocs test.
%    01/24/18  jnm  Formatting update to match Wiki.

if notDefined('dataType'), dataType = 'photons'; end
if notDefined('roiLocs'), error('roi locs required'); end

wave = oiGet(oi, 'wave');
irradiance = vcGetROIData(oi, roiLocs, dataType);
irradiance = mean(irradiance);

if length(wave) == 1
    % For a monochrome image, a plot doesn't make any sense. So, we just
    % put up a box describing the mean irradiance.
    switch dataType
        case 'photons'
            str = sprintf('Irradiance: %.3e (q/s/m^2/nm)  at %.0f nm', ...
                irradiance, wave);
        case 'energy'
            str = sprintf('Irradiance: %.3e (Watts/m^2/nm) at %.0f nm', ...
                irradiance, wave);
        otherwise
            error('Unknown data type.');
    end
    msgbox(str);
else
    % Attach data to the figure itself
    udata.x = wave;
    udata.y = irradiance;
    udata.roiLocs = roiLocs;

    plot(wave, irradiance);
    set(gca, 'ylim', ...
        [.95 * min(irradiance(:)), 1.03 * max(irradiance(:))]);
    xlabel('Wavelength (nm)');
    grid on;

    switch lower(dataType)
        case 'photons'
            ylabel('Irradiance (q/s/m^2/nm)');
        case 'energy'
            ylabel('Irradiance (Watts/m^2/nm)');
        otherwise
            disp('Unknown data type')
    end
end

end

function uData = oiPlotOTF(oi, pType, varargin)
% Plot OTF functions associated with the optics in an optical image
%
% Syntax:
%   oiPlotOTF([oi], [pType])
%
% Description:
%    Plot OTF functions associated with the optics in an optical image.
%
% Inputs:
%    oi    - (Optional) The optical image. Default is  to query vcGetObject
%    pType - (Optional) The Plot Type. Default 'otf550' Options include:
%       {'otf'}           - Optical transfer function, units are lines/mm
%       {'otf550'}        - OTF at 550 nm
%       {'psf'}           - Point spread function at selected wavelength
%       {'psf550'}        - PSF at 550nm spatial units are microns
%       {'lswavelength'}  - Line spread function at all wavelengths.
%                           Peak spatial frequency can be set for the OTF
%                           (default: 3 * incoherent cutoff). Number of
%                           spatial samples to plot in the line spread can
%                           be set (default: 40).
%       {'otfwavelength'} - One dimensional cut through the OTF at a all
%                           wavelengths. Units are cycles/mm
%
% Outputs:
%    uData - The user data structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%	 * [Note: JNM - When 2015b cycles out, replace strfind with contains]
%    * [Note: XXX - TODO: Implement raytrace
%    * [Note: JNM - TODO: Someone needs to check over the "Note: XXX - "
%      notes below, so we can determine which are safe to remove. There are
%      several such notes.]
%    * [Note: XXX: Determine how to better select the number of samples for
%      the spatial frequency. Currently 100 samples, the number of which is
%      arbitrarily chosen.]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    xx/xx/12       Moved into oiPlot June, 2012.
%    12/12/17  jnm  Formatting
%    01/24/18  jnm  Formatting update to match Wiki.

if notDefined('oi'), oi = vcGetObject('oi'); end
if notDefined('pType'), pType = 'otf550'; end

wavelength = oiGet(oi, 'wavelength');
optics = oiGet(oi, 'optics');

% This catches the case in which the oi has not yet been defined, but the
% optics have.
if isempty(wavelength)
    oi = initDefaultSpectrum(oi, 'hyperspectral');
    optics = initDefaultSpectrum(optics, 'hyperspectral');
    wavelength = opticsGet(optics, 'wavelength');
end

nWave = oiGet(oi, 'nwave');
pType = ieParamFormat(pType);

switch lower(pType)
    case {'otf', 'otf550'}
        % oiPlotOTF(oi, 'otf', thisWave, nSamp);
        % OTF at a selected wavelength.
        units = 'mm';  % Units are cycles/mm
        if strfind(pType, '550') %#ok<*STRIFCND>
            thisWave = 550;
        elseif length(varargin) >= 1
            thisWave = varargin{1};
        else
            thisWave = ieReadNumber('Select OTF wavelength (nm)', 550, ...
                '%.0f');
        end

        if length(varargin) >= 2
            nSamp = varargin{2};
        else
            nSamp = 40;
        end

        % Retrieve OTF data (which might be complex) from the optics
        opticsModel = opticsGet(optics, 'opticsModel');
        switch lower(opticsModel)
            case {'dlmtf', 'diffractionlimited'}
                % Compute the otf data

                % Specify frequency support and compute the dl MTF. Note
                % that DC is in the center of the support vectors, while
                % OTF returned by dlMTF has DC at the (1,1) upper left
                % position.
                %
                % Over sample the frequency support just to make it big
                freqOverSample = 2;
                fSupport = opticsGet(optics, 'dl fsupport matrix', ...
                    thisWave, units, nSamp);

                fSupport = fSupport * freqOverSample;
                otf = dlMTF(oi, fSupport, thisWave, units);

                % DC is at (1, 1) in the returned OTF; we plot with DC in
                % the center.
                otf = fftshift(otf);
                figTitle = sprintf('DL OTF at %.0f', thisWave);

            case {'shiftinvariant'}
                % Get OTF from optics structure.  Error if it wasn't there.
                otf = opticsGet(optics, 'otf data', thisWave);
                if isempty(otf), error('No OTF data'); end

                % Units are cycles/mm of optics support.  Note that DC is
                % in the center of the support vectors, while OTF returned
                % below has DC at the (1,1) upper left position.
                s = opticsGet(optics, 'otfSupport');
                fSupport(:, :, 1) = s{1};
                fSupport(:, :, 2) = s{2};

                % Transform so DC is in center for plotting.
                otf = fftshift(otf);
                figTitle = sprintf('abs(OTF) at %.0f nm', thisWave);

            case {'raytrace'}
                error('Ray trace plot: Not yet implemented');

            otherwise
                error('Unknown optics model: %s\n', opticsModel);
        end

        % Select the support and plot the mesh
        sz = selectPlotSupport(otf, 0.01);
        x = getMiddleMatrix(fSupport(:, :, 1), sz);
        y = getMiddleMatrix(fSupport(:, :, 2), sz);
        otf = getMiddleMatrix(otf, sz);
        if isreal(otf(:))
            mesh(x, y, otf);
        else
            disp('Complex otf values'), mesh(x, y, abs(otf));
        end

        % Label axes and store data
        xlabel('cyc/mm');
        ylabel('cyc/mm');
        zlabel('amplitude');
        title(figTitle);
        uData.otf = otf;
        uData.fSupport = fSupport;

    case {'psf', 'psf550'}
        % Spatial scale is microns.
        units = 'um';
        nSamp = 100;
        freqOverSample = 4;
        if strfind(pType, '550')
            thisWave = 550;
        elseif length(varargin) >= 1
            thisWave = varargin{1};
        else
            thisWave = ieReadNumber('Select PSF wavelength (nm)', 550, ...
                '%.0f');
        end

        opticsModel = opticsGet(optics, 'model');
        switch lower(opticsModel)
            case {'diffractionlimited'}
                % The opticsGet() for diffraction limited should be
                % adjusted so that this code becomes shorter.

                psf = opticsGet(optics, 'psf data', ...
                    thisWave, units, nSamp, freqOverSample);

                fSupport = opticsGet(optics, 'dlFSupport matrix', ...
                    thisWave, units, nSamp);
                fSupport = fSupport * freqOverSample;

                % Put samples symmetric around 0
                % Make them spaced properly
                deltaSpace = 1 / (2 * max(fSupport(:)));
                samp = (-nSamp:(nSamp - 1));
                [X, Y] = meshgrid(samp, samp);
                sSupport(:, :, 1) = X * deltaSpace;
                sSupport(:, :, 2) = Y * deltaSpace;

                % Plot a black circle at the first zero of the Airy disk.
                % Convert the wavelength to meters. Multiply by the Airy
                % disk formula to determine the first zero crossing. Scale
                % units to microns, because that is how we plot the data.
                fNumber = opticsGet(optics, 'fNumber');
                radius = (2.44 * fNumber * thisWave * 10 ^ -9) / 2 ...
                    * ieUnitScaleFactor('um');
                nSamp = 200;
                [~, ptsXY] = ieShape('circle', 'nSamp', nSamp, ...
                    'radius', radius);
                adX = ptsXY(:, 1);
                adY = ptsXY(:, 2);
                adZ = zeros(size(ptsXY(:, 1)));
            case {'shiftinvariant'}
                psf = opticsGet(optics, 'psfdata', thisWave);
                tmp = opticsGet(optics, 'psfSupport', units);
                sSupport(:, :, 1) = tmp{1};
                sSupport(:, :, 2) = tmp{2};

            case {'raytrace'}
                % opticsGet(optics, 'rtPSFdata') should be cleaned up for
                % this call. Spatial support, frequency support, all of
                % that should be in there.
                error('Not yet implemented');

            otherwise
                error('Unknown otf function: %s\n', opticsModel);
        end

        % Plot it and if DL, then add the Airy disk
        mesh(sSupport(:, :, 1), sSupport(:, :, 2), psf);
        if strcmpi(opticsModel, 'diffractionlimited')
            hold on;
            plot3(adX, adY, adZ, 'k.');
            hold off;
        end

        % Label, store data
        xlabel('Position (um)');
        ylabel('Position (um)');
        zlabel('Irradiance (relative)');
        title(sprintf('Point spread (%.0f nm)', thisWave));
        uData.x = sSupport(:, :, 1);
        uData.y = sSupport(:, :, 2);
        uData.psf = psf;

    case {'lswavelength'}
        % Line spread function at all wavelengths
        units = 'um';
        wavelength = oiGet(oi, 'wavelength');
        nWave = oiGet(oi, 'nwave');
        model = opticsGet(optics, 'model');

        % Choose the peak frequency for the OTF. If none is passed in, we
        % use the incoherent cutoff frequency.
        if length(varargin) >= 1 && ~isempty(varargin{1})
            peakF = varargin{1};
        else
            switch lower(model)
                case 'diffractionlimited'
                    inCutoff = opticsGet(optics, 'inCutoff', units);
                    peakF = 3 * max(inCutoff);
                case 'shiftinvariant'
                    fx = opticsGet(optics, 'otffx', 'um');
                    peakF = max(abs(fx(:)));
                otherwise
                    error('LS not implemented for %s model', model);
            end
        end
        if length(varargin) >= 2
            spaceSamp = varargin{2};
        else
            spaceSamp = 40;
        end

        % [Note: XXX - The spaceSamp and nSamp parameters are not clearly
        % enough defined. The reason we care is because the code is broken
        % when spaceSamp is not 40. The problem appears to be that lsfWave
        % computed below might have only 60 samples and we might ask for,
        % say 120. So, we should at least check.]

        % The incoherent cutoff frequency has units of cycles/micron
        % So, 1/inCutoff has units of microns/Nyquist
        % The maximum frequency is at the Nyquist, and there are two
        % samples at the Nqyuist. So the sample spacing is half the peakF
        deltaSpace = 1 / (2 * peakF);

        % Make the spatial frequency samples used to compute the OTF. These
        % run from [-peakF, +peakF]. We make 100 samples, which is pretty
        % arbitrary. Not sure how to choose this better. Should be using
        % unitFrequencyList() here.
        nSamp = 100;
        fSamp = (-nSamp:(nSamp - 1)) / nSamp;
        [fX, fY] = meshgrid(fSamp, fSamp);
        fSupport(:, :, 1) = fX * peakF;
        fSupport(:, :, 2) = fY * peakF;

        opticsModel = opticsGet(optics, 'opticsModel');
        switch lower(opticsModel)
            case {'dlmtf', 'diffractionlimited'}
                otf = dlMTF(oi, fSupport, wavelength, units);
            case {'shiftinvariant'}
                sz = opticsGet(optics, 'otf size');
                otf = zeros(sz(1), sz(2), nWave);
                for ii = 1:nWave
                    otf(:, :, ii) = opticsGet(optics, 'otfdata', ...
                        wavelength(ii));
                end
            otherwise
                error('LSWavelength1D not implemented for model: %s\n', ...
                    opticsModel);
        end

        % Create the line spread for a horizontal line. We use the first
        % row of the otf to estimate the line spread. This only works if
        % the OTF is circularly symmetric; if it is not, there isn't really
        % a single line spread.

        % [Note: XXX - We should figure out the right value of spaceSamp
        % In some cases spaceSamp is bigger than the otf. So we need to
        % cut it back]
        spaceSamp = min(spaceSamp, size(otf, 2) - 1);
        % lsWave = zeros(spaceSamp, nWave);
        for ii = 1:nWave
            % The central line in the otf is the first line
            tmp = otf(1, :, ii);  % figure; imagesc(abs(otf(:, :, ii)))

            % We invert the OTF along that line to get an LSF. We apply
            % the fftshift because we want the lsf to be centered.
            lsf = fftshift(ifft(tmp));   % figure; plot(abs(lsf))

            % Pull out samples from the middle because otherwise the image
            % can be hard to see.
            lsWave(:, ii) = getMiddleMatrix(lsf, spaceSamp);
        end

        % Choose the x coords that match the line spread spatial samples.
        X = (-nSamp:(nSamp - 1)) * deltaSpace;
        X = getMiddleMatrix(X, spaceSamp);

        % Show it
        if nWave > 1
            mesh(X, wavelength, lsWave');
            xlabel('Position (um)');
            ylabel('Wavelength (nm)');
            zlabel('Intensity (rel.)');
        else
            plot(X, lsWave(:));
            xlabel('Wavelength (nm)');
            ylabel('Intensity (rel.)');
        end
        view(30, 20);

        % Store the results in the figure.
        uData.x = X;
        uData.wavelength = wavelength;
        uData.lsWave = lsWave';

    case {'otfwavelength'}
        % Plot a line through the center of the OTF as a function of
        % wavelength
        opticsModel = opticsGet(optics, 'opticsModel');
        units = 'um';

        % [Note: BW - We get the OTF slightly differently for the different
        % models. If we rewrote opticsGet to check for the optics model, we
        % could do things a little more simply here. Maybe we should put
        % this code into opticsGet.]
        switch lower(opticsModel)
            case {'dlmtf', 'diffractionlimited'}
                % Make the spatial frequency samples used to compute the
                % OTF. These run from [-peakF, +peakF]. We make 100
                % samples, which is pretty arbitrary. Not sure how to
                % choose this better. Should be using unitFrequencyList()
                % here.
                %
                % [Note: BW - Determine how to better select the number of
                % samples for the spatial frequency.]
                if length(varargin) >= 1
                    peakF = varargin{1};
                else
                    inCutoff = opticsGet(optics, 'inCutoff', units);
                    peakF = 3 * max(inCutoff);
                end
                nSamp = 100;
                fSamp = (-nSamp:(nSamp - 1)) / nSamp;
                [fX, fY] = meshgrid(fSamp, fSamp);
                fSupport(:, :, 1) = fX * peakF;
                fSupport(:, :, 2) = fY * peakF;

                % Here the OTF has DC at (1,1).  Will deal with that below.
                otf = dlMTF(oi, fSupport, wavelength, units);

            case {'shiftinvariant'}
                % Data are stored in OTF slot, with DC at (1,1).  Will deal
                % with that below.
                s = opticsGet(optics, 'otfSupport');
                fSupport(:, :, 1) = s{1};
                fSupport(:, :, 2) = s{2};
                otf = zeros(length(s{1}), length(s{2}), nWave);
                for ii = 1:nWave
                    otf(:, :, ii) = abs(opticsGet(optics, 'otfdata', ...
                        wavelength(ii)));
                end
            otherwise
                error('OTF 1D plot not implemented for: %s\n', ...
                    opticsModel);
        end

        fx = fSupport(1, :, 1);
        otfWave = zeros(length(fx), nWave);

        % Convert from DC at (1,1) to DC at center, using fftshift.
        for ii = 1:nWave, otfWave(:, ii) = fftshift(otf(1, :, ii)); end

        mesh(fx, wavelength, otfWave');
        view(30, 20);
        xlabel('cycles/mm');
        ylabel('Wavelength (nm)');
        zlabel('abs(OTF)');

        % Store the data in the figure.
        uData.otf = otfWave;
        uData.fSupport = fx;
        uData.wavelength = wavelength;

    otherwise
        error('Unknown oiPlotOTF data type.');
end

end

function uData = oiPlotIlluminanceMesh(oi, yScale)
% Plot optical image illuminance (lux) as a mesh
%
% Syntax:
%   oiPlotIlluminanceMesh(oi, yScale)
%
% Description:
%    The default scaling of the lux axis is logarithmic. Set yScale to
%    'linear' for a linear scale.
%
%    If roiFlag is set to true (1), the user selects a region of the image
%    from the optical image window. Otherwise the entire image is plotted.
%
%    In the GraphWin, set Tools | Move Camera to rotate the view.
%
% Inputs:
%    oi     - The optical image.
%    yScale - The scaling method for the y-axis. Default is 'log' Options
%             are 'log' and 'linear'.
%
% Outputs:
%    uData  - The user data structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - the variable roiFlag is not included in the header, but
%      there is a query for whether or not it was defined at the beginning
%      of the function? It is also not used, so I removed the line.]
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('oi'), error('OI required.'); end
if notDefined('yScale'), yScale = 'log'; end

illum = oiGet(oi, 'illuminance');

spacing = oiGet(oi, 'sample spacing', 'um');
sz = size(illum);
r = (1:sz(1)) * spacing(1);
c = (1:sz(2)) * spacing(2);
switch yScale
    case 'log'
        uData.data = fliplr(log10(illum));
        mesh(c, r, uData.data);
        zlabel('Lux (Log 10)')
    case 'linear'
        uData.data = fliplr(illum);
        mesh(c, r, uData.data);
        zlabel('Lux')
    otherwise
        error('unknown yScale.');
end

uData.c = c;
uData.r = r;

xlabel('um');
ylabel('um');
title('Illuminance');

end

function uData = oiPlotCIE(oi, dataType, roiLocs)
% Plot CIE data from optical image.
%
% Syntax:
%   uData = oiPlotCIE(oi, dataType, roiLocs)
%
% Description:
%    Could be moved into the case statements of the mother ship.
%
%    Graph  optical image properties (Luminance, chromaticity coordinates)
%    from an ROI. The user is prompted to select the ROI in the OI window.
%
%    The plotted values can be obtained from the GraphWin using
%
%       udata = get(gcf, 'userdata');
%
% Inputs:
%    oi       - The optical image
%    dataType - The type of oi plot. Options are: chromacity, illuminance.
%    roiLocs  - The region of interest locations
%
% Outputs:
%    uData    - The user data structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - There are no default options for any of the arguments,
%      or checks that they are populated.]

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/12/17  jnm  Formatting

% Examples:
%{
    % [Note: JNM - Example will not work as this is a nested function and
    % not publicly accessible. Calling oiPlot with a plot type of
	% chromaticityroi, however, will end up using this function.]
    oi = vcGetObject('oi');
    vcNewGraphWin;
    udata = oiPlotCIE(oi, 'chromaticity')
    oiPlotCIE(oi, 'illuminance', roiLocs);
    oiPlotCIE(oi, 'chromaticity', roiLocs);
%}

switch lower(dataType)
    case {'chromaticity'}
        photons = vcGetROIData(oi, roiLocs, 'photons');
        wave = oiGet(oi, 'wave');
        XYZ = ieXYZFromPhotons(photons, wave);
        data = chromaticity(XYZ);
        uData.x = data(:, 1);
        uData.y = data(:, 2);
        val = mean(XYZ);
        valxy = mean(data);

        chromaticityPlot(data, [], [], 0);
        title('roiLocs-chromaticity (CIE 1931)');

        txt = sprintf('Means\n');
        tmp = sprintf('X = %.02f\nY = %.02f\nZ = %.02f\n', ...
            val(1), val(2), val(3));
        txt = addText(txt, tmp);
        tmp = sprintf('x = %0.02f\ny = %0.02f\n', valxy(1), valxy(2));
        txt = addText(txt, tmp);
        text(0.8, 0.65, txt);
        axis equal
        hold off

    case {'illuminance'}
        data = vcGetROIData(oi, roiLocs, 'illuminance');
        hist(data(:));
        uData.illum = data;
        xlabel('Iluminance (lux)');
        ylabel('Count');
        title('Iluminance histogram');

    otherwise
        error('Unknown oi plot data type %s\n', dataType);
end

uData.roiLocs = roiLocs;
oName = oiGet(oi, 'name');
set(gcf, 'Name', sprintf('ISET-OI: %s', oName));

end


function sz = selectPlotSupport(data, prct)
% Used with getMiddleMatrix to pull out the 'interesting' center of a plot
%
% Syntax:
%   sz = selectPlotSupport(data, prct)
%
% Description:
%    Sometimes we have a large surface to plot but the interesting
%    part is near the middle of the data set. Rather than plotting the
%    entire surf or mesh(data) we pull out a central region. This  is
%    the method for choosing the  size of the data we pull out. This
%    method is used in conjunction with getMiddleMatrix.
%
% Inputs:
%    data - The data set to plot
%    prct - (Optional) What percentage of the middle you wish to display,
%           in decimal form. Default is 0.01 (1%)
%
% Outputs:
%    sz   - The center value(s) of data of size prct
%
% Optional key/value pairs:
%    None.
%
% Notes:
%	 * [Note: XXX - What if data are a vector? Can we adjust this routine
%	   to make it work?]
%
%  See Also:
%    meshPlot, plotOTF

% Examples:
%{
 g = fspecial('gaussian',256,30);
 g = ieScale(g,0,1);
 r = selectPlotSupport(g,0.1);
 vcNewGraphWin; imagesc(g); axis image
 ieShape('circle','center',[128,128],'radius',r,'color','white');
 colorbar;
%}

%%

if notDefined('prct'), prct = 0.01; end
minRadius = 25;

r  = size(data, 1);
mx = max(data(:));
centerRow = round(r / 2);

%% Find the location in the center row at the edge of the prct max
l = (data(centerRow, :) < prct * mx);

if max(l) == 0, radius = centerRow - 1;
else
    % The returned value is at least 25.
    [~, idx] = max(data(centerRow, l));
    radius = centerRow - idx;
    if radius < 25
        warning('Selected radius is very small. Increasing to %d',minRadius);
    end
end
sz = radius;
end