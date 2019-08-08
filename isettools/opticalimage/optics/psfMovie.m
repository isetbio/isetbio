function psfMovie(oi, figNum, delay)
% Show a movie of the pointspread functions
%
% Syntax:
%   psfMovie([oi], [figNum], [delay])
%
% Description:
%    The movies differ slightly for the shift-invariant and ray trace
%    methods. The shift-invariant doesn't depend on field height; the ray
%    trace does. So we show the full set differently.
%
%    For the shift invariant we just show a movie over wavelength
%
%    For the ray trace
%    The image height increase along the x-axis. The data are displayed for
%    one wavelength first, and then the next wavelength. The horizontal
%    axis indicates the image height in microns.
%
% Inputs:
%    oi     - (Optional) Struct. An optical image structure. Default uses
%             vcGetObject to retrieve an existing instance.
%    figNum - (Optional) Numeric. Figure number. Default to pulling an
%             existing figure using vcSelectFigure.
%    delay  - (Optional) Numeric. Delay in seconds. Default 0.2.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval, LLC, 2005
%    03/13/18  jnm  Formatting
%    04/07/18  dhb  Add ETTBSkip for broken examples.
%    06/21/18  dhb  Someone added more broken examples, or they have
%                   broken since the last time I was here. Now they
%                   are skipped.
%    06/26/19  JNM  Documentation update. Notes on failing examples inline.

% Examples:
%{
    % ETTBSkip.  Example broken.
    %
    % vcGetObject retursn empty matrix
    oi = vcGetObject('oi');
    psfMovie(oi, 1)
%}
%{
    % ETTBSkip.  Example broken.
    %
    % This example dies on the call to vcImportObject
    vcImportObject('OPTICS');
    psfMovie([], 1);
    optics = vcGetObject('optics');
    optics = rtPSFEdit(optics, 0, 1, 2);

    % Returns error: "Error using chdir
    % Cannot CD to session (Name is nonexistent or not a directory)."
    % [Note: JNM - rtPSFEdit does not exist?]
%}
%{
    % ETTBSkip.  Example broken.
    %
    % Crashes on a call to opticsGet.
    oi = oiCreate('shift invariant');
    psfMovie(oi);
    % [Note: JNM - unknown optics parameter 'shiftinvariantpsfdata'. L126.]
%}
%{
    % ETTBSkip.  Example broken.
    %
    % Crashes on a call to opticsGet.
    oi = oiCreate('diffraction limited');
    psfMovie(oi);
    % [Note: JNM - same error as shift invariant (psf data call) but L111.]
%}

%%
if notDefined('oi'), oi = vcGetObject('oi'); end
if notDefined('figNum'), figNum = vcSelectFigure('GRAPHWIN'); end
if notDefined('delay'), delay = 0.2; end

figure(figNum)
set(figNum, 'name', 'PSF Movie');

optics = oiGet(oi, 'optics');
opticsModel = oiGet(oi, 'optics model');

switch lower(opticsModel)
    case 'diffractionlimited'
        % This needs to be checked. The psf doesn't seem to be changing
        % correctly with wavelength.
        %
        % The OTF is computed on the fly for the diffraction limited case.
        %
        nSamp = 16;
        freqOverSample = 4;
        fSupport = oiGet(oi, 'fsupport', 'um') * freqOverSample;
        samp = (-nSamp:(nSamp - 1));
        % [X, Y] = meshgrid(samp, samp);
        % deltaSpace = 1 / (2 * max(fSupport(:)));
        % sSupport(:, :, 2) = Y * deltaSpace;
        % sSupport(:, :, 1) = X * deltaSpace;

        deltaSpace = 1 / (2 * max(fSupport(:)));
        x = samp * deltaSpace;
        y = x;
        wave = opticsGet(optics, 'wave');

        vcNewGraphWin;
        for ii=1:length(wave)
            psf = opticsGet(optics, 'diffraction limited psf data', ...
                wave(ii), 'um', nSamp, freqOverSample);
            %{
            imagesc(y, x, psf(:, :));
            xlabel('Position (um)');
            ylabel('Position (um)');
            grid on;
            axis image;
            title(sprintf('Wave %.0f nm', wave(ii)));
            pause(delay);
            %}
        end

    case 'shiftinvariant'
        % We get the psf data all at once in this case
        psf = opticsGet(optics, 'shift invariant psf data');
        support = opticsGet(optics, 'psf support', 'um');
        y = support{1}(:); x = support{2}(:);
        wave = opticsGet(optics, 'wavelength');
        w = size(psf, 3);

        for ii=1:w
            imagesc(y, x, psf(:, :, ii));
            xlabel('Position (um)');
            ylabel('Position (um)');
            grid on;
            axis image
            title(sprintf('Wave %.0f nm', wave(ii)));
            pause(delay);
        end
    %{
    case 'raytrace'
        name = opticsGet(optics, 'rtname');
        figNum = vcNewGraphWin;
        set(figNum, 'name', sprintf('%s: PSF movie', name));
        colormap(gray(256));

        wave = opticsGet(optics, 'rt psf wavelength');
        imgHgt = opticsGet(optics, 'rt psf field height', 'um');
        psf = opticsGet(optics, 'rt psf data');
        c = opticsGet(optics, 'rt psf support col', 'um');
        r = opticsGet(optics, 'rt psf support col', 'um');

        % Should we plot them on a single image and move them, or centered
        % like this?
        gColor = [.5 .5 0];
        for jj = 1:length(wave)
            for ii = 1:length(imgHgt)
                imagesc(r + imgHgt(ii), c + imgHgt(ii), ...
                    squeeze(psf(:, :, ii, jj)));
                set(gca, 'yticklabel', []);
                xlabel('Position (um)');
                set(gca, 'xcolor', gColor, 'ycolor', gColor);
                grid on;
                axis image
                title(sprintf('Wave %.0f nm\nField height %.2f um', ...
                    wave(jj), imgHgt(ii)));
                pause(delay);
            end
        end
    %}
    otherwise
        error('Unknown model %s\n', opticsModel);
end

end