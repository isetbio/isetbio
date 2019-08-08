function optics = siSynthetic(psfType, oi, varargin)
% Create synthetic shift-invariant optics and insert into optics structure
%
% Syntax:
%   optics = siSynthetic(psfType, oi, [varargin])
%
% Description:
%    This code was used for testing the shift-invariant optics.  We build
%    on this to let the user create a custom shift-invariant optics.
%
%    By default, the optics (custom) fields are filled in using simple
%    values.
%
%    Must keep spatial sampling parameters pretty consistent across usages.
%
%    This function and its cousins may not be around for long in ISETBIO.
%
%    Examples are included within the code. To access, type 'edit
%    siSynthetic.m' into the Command Window.
%
% Inputs:
%    psfType  - String. The PSF Type. Options are:
%         'gaussian' --  bivariate normals.
%         'custom'   --  read a file with variables explained below
%    oi       - Struct. An optical image. It must be shift invariant type
%    varargin - (Optional) VARIES. Additional arguments for the function.
%               Arguments depend on the psf type. See some examples below.
%       for gaussian:
%           waveSpread: size of the PSF spread at each of the wavelength
%                       for gaussian this is in microns (um)
%           xyRatio:    Ratio of spread in x and y directions
%           filename:   Output file name for the optics
%
%       for custom:
%           inData  - filename or struct with psf, umPerSamp, and wave data
%           outFile - Optional
%
% Outputs:
%    optics   - Struct. The created optics structure.
%
% Optional key/value pairs:
%   None.
%
% See Also:
%   s_SIExamples, ieSaveSIOpticsFile, t_codeFFTinMatlab
%

% History:
%    xx/xx/05        Copyright ImagEval Consultants, LLC, 2005.
%    12/08/17  dhb  Take number of otf samples from oi, not hard code @ 128
%              dhb  Take mm/[psf sample] from oi, not hard code at 0.25e-3.
%              BW   We might just eliminate a lot of this set of SI methods
%              dhb  Created working example for 'custom'
%    01/22/18  dhb  Examples run in clean workspace.
%    03/15/18  jnm  Formatting
%    06/28/19  JNM  Documentation update

% Examples:
%{
    oi = oiCreate('shift invariant');
    oi = initDefaultSpectrum(oi, 'spectral');
    wave = oiGet(oi, 'wave');
    psfType = 'gaussian';
    waveSpread = wave / wave(1);

    % Circularly symmetric Gaussian
    xyRatio = ones(1, length(wave));
    optics = siSynthetic(psfType, oi, waveSpread, xyRatio);
    oi = oiSet(oi, 'optics', optics);
    oiPlot(oi, 'otf 550');
%}
%{
    oi = oiCreate('shift invariant');
    oi = initDefaultSpectrum(oi, 'spectral');
    wave = oiGet(oi, 'wave');
    psfType = 'gaussian';
    waveSpread = wave / wave(1);

    % Make one with an asymmetric Gaussian
    xyRatio = 2 * ones(1, length(wave));
    optics =  siSynthetic(psfType, oi, waveSpread, xyRatio);
    oi = oiSet(oi, 'optics', optics);
    oiPlot(oi, 'otf 550');
%}
%{
    % Create oi and get info about OTF/PSF dimensions and support
    oi = oiCreate;
    optics = oiGet(oi, 'optics');
    dx(1:2) = opticsGet(optics, 'psf spacing', 'um');
    fx = opticsGet(optics, 'otf fx');
    nSamples = length(fx);

    % Create an SI data file and then use it. We create a delta function
    % PSF.  Need to match what is in oi for siSynthetic to work.
    psf = zeros(nSamples, nSamples, 31);
    psf(floor(nSamples / 2) + 1, floor(nSamples / 2) + 1, :) = 1;
    wave = 400:10:700;
    umPerSamp = dx;
    fname = tempname;
    ieSaveSIDataFile(psf, wave, umPerSamp, fname);

    % Use siSynthetic to get optics and put those into the oi.
    optics = siSynthetic('custom', oi, fname);
    oi = oiSet(oi, 'optics', optics);
    delete([fname '.mat']);

    % Since PSF was delta, should get all 1's in OTF
    oiPlot(oi, 'otf 550');
%}

%% Parameter initializiation
if notDefined('psfType'), psfType = 'gaussian'; end
if notDefined('oi'), oi = vcGetObject('oi'); end

if ~strcmp(oiGet(oi, 'optics model'), 'shiftinvariant')
    error('oi model must be "shift invariant"');
end

inFile = [];
outFile = [];

% Wavelength samples
wave = oiGet(oi, 'wave');
nWave = length(wave);

%% Make array for new OTF
optics = oiGet(oi, 'optics');
otfSamples = opticsGet(optics, 'otf size');
if (otfSamples(1) ~= otfSamples(2))
    error('OTF must be on square support');
end
nSamples = otfSamples(2);
OTF = zeros(nSamples, nSamples, nWave);

% Get spacing of psf in mm of retina
dx(1:2) = opticsGet(optics, 'psf spacing', 'mm');

%% Create PSF and OTF
switch lower(psfType)
    case 'gaussian'
        % Create a Gaussian set of PSFs.
        if length(varargin) < 2
            error('Wavespread and xyRatio required');
        end
        xSpread = varargin{1};    % Spread is in units of um here
        xyRatio = varargin{2};
        ySpread = xSpread(:) .* xyRatio(:);
        if length(varargin) == 3, outFile = varargin{3};
        else, outFile = [];
        end

        % Convert spread from microns to millimeters because OTF data are
        % stored in line pairs per mm
        xSpread = xSpread / 1000;
        ySpread = ySpread / 1000;

        for jj = 1:nWave
            % We convert from spread in mm to spread in samples for the
            % biNormal calculation.
            psf = biNormal(xSpread(jj) / dx(2), ySpread(jj) / dx(1), ...
                0, nSamples);
            psf = psf / sum(psf(:));

            % Use PsfToOtf to make the change, and then put center in upper
            % right to match isetbio conventions.  Commented out below is
            % the older code, which may or may not do the same thing
            [~,~,centeredOTF] = PsfToOtf([],[],psf);
            OTF(:,:,jj) = ifftshift(centeredOTF);
            % psf = fftshift(psf);
            % OTF(:, :, jj) = fft2(psf);
        end

    case 'custom'
        %% Get PSF data
        if isempty(varargin)
            % Find a file by asking user
            inFile = ...
                vcSelectDataFile('stayPut', 'r', 'mat', ...
                'Select custom SI optics');
            if isempty(inFile)
                disp('User canceled');
                optics = [];
                return;
            end
        elseif ischar(varargin{1})
            % This is a file name.  Load it and get parameters
            tmp = load(varargin{1});
            if ~isfield(tmp, 'psf'), error('Missing psf variable');
            else, psfIn = tmp.psf; end
            if ~isfield(tmp, 'wave'), error('Missing wave variable');
            else, wave = tmp.wave; end
            if ~isfield(tmp, 'umPerSamp'), error('Missing wave variable');
            else, mmPerSamp = (tmp.umPerSamp) / 1000; end
        elseif isstruct(varargin{1})
            % The data were sent in as a iset shift-invariat PSF struct
            psfIn = varargin{1}.psf;
            wave = varargin{1}.wave;
            mmPerSamp = (varargin{1}.umPerSamp) / 1000;
        end
        if length(varargin) > 1, outFile = varargin{2}; end

        % Check the parameters for consistency
        [m, n, nWave] = size(psfIn);
        if length(wave) ~= nWave
            error('Mis-match between wavelength and psf');
        end
        if m ~= nSamples || n ~= nSamples
            error('Need input and output number of samples to match');
        end
        if (max(abs(mmPerSamp(2) - dx(2))) > 1e-10 || ...
                max(abs(mmPerSamp(1) - dx(1))) > 1e-10)
            error('Cannot yet change psf sampling here.')
        end

        % OTF computation
        %
        % This is the sampling grid of the psfIn.
        % Units at this point are in mm. The psf gets interpolated
        % to the desired sampling size and then converted to an OTF.
        x = (1:n) * mmPerSamp(2);
        x = x - mean(x(:));
        y = (1:m) * mmPerSamp(1);
        y = y - mean(y(:));
        [xInGrid, yInGrid] = meshgrid(x, y);

        xOut = (1:nSamples) * dx(2);
        xOut = xOut - mean(xOut(:));
        yOut = (1:nSamples) * dx(1);
        yOut = yOut - mean(yOut(:));
        [xOutGrid, yOutGrid] = meshgrid(xOut, yOut);

        for ii = 1:nWave
            psf = interp2(xInGrid, yInGrid, psfIn(:, :, ii), xOutGrid, ...
                yOutGrid, 'linear', 0);
            psf = psf / sum(psf(:));

            % Use PsfToOtf to make the change, and then put center in upper
            % right to match isetbio conventions.  Commented out below is
            % the older code, which may or may not do the same thing
            [~, ~, centeredOTF] = PsfToOtf([], [], psf);
            OTF(:, :, ii) = ifftshift(centeredOTF);
            %psf = fftshift(psf);
            %OTF(:, :, ii) = fft2(psf);
        end

    otherwise
        error('Unspecified PSF format');
end

%% Find OI sample spacing.  The OTF line spacing is managed in lines/mm
nyquistF = 1 ./ (2 * dx);   % Line pairs (cycles) per mm
fx = unitFrequencyList(nSamples) * nyquistF(2);
fy = unitFrequencyList(nSamples) * nyquistF(1);

% [FY, FX] = meshgrid(fy, fx);
% vcNewGraphWin;
% mesh(FY, FX, fftshift(abs(OTF(:, :, 2))))

%% Create and save the optics
optics = opticsCreate;
optics = opticsSet(optics, 'model', 'shiftinvariant');
optics = opticsSet(optics, 'name', inFile);
optics = opticsSet(optics, 'otffunction', 'custom');
optics = opticsSet(optics, 'otfData', OTF);
optics = opticsSet(optics, 'otffx', fx);
optics = opticsSet(optics, 'otffy', fy);
optics = opticsSet(optics, 'otfwave', wave);
optics.lens.density = 0;

if isempty(outFile), return; else, vcSaveObject(optics, outFile); end

end