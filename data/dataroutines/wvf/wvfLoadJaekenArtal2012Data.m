function [wvf, oi] = wvfLoadJaekenArtal2012Data(varargin)

%
% Syntax:
%    zCoefs = wvfLoadJaekenArtal2012Data(jIndex, whichEye, eccentricity)
%
% Description:
%    Convert a dataset of individual patient Zernike coefficients
%    to a mean OTF for eccentricities -40 to 40 degrees on the horizontal axis
%    for a given eye. Pupil size was 4 mm.
%
%    Data are published in Jaeken, B. & Artal, P. (2012) Optical Quality of
%    Emmetropic and Myopic Eyes in the Periphery Measured with High-Angular
%    Resolution. Investigative Ophthalmology & Visual Science, June 2012,
%    Vol. 53, No. 7
%
%    For information on Zernike coefficient and their names:
%	 http://www.telescope-optics.net/monochromatic_eye_aberrations.htm
%
%    Table of names
%    =================================
%      j   name
%
%      0  'piston'
%      1  'vertical_tilt'
%      2  'horizontal_tilt'
%      3  'oblique_astigmatism'
%      4  'defocus'
%      5  'vertical_astigmatism'
%      6  'vertical_trefoil'
%      7  'vertical_coma'
%      8  'horizontal_coma'
%      9  'oblique_trefoil'
%      10 'oblique_quadrafoil'
%      11 'oblique_secondary_astigmatism'
%      12 'primary_spherical', 'spherical'
%      13 'vertical_secondary_astigmatism'
%      14 'vertical_quadrafoil'
%
% Inputs:
%    jIndex       - List of OSA J values
%    whichEye     - string defining whether to use left or right eye
%    eccentricity - integer between -40 and 40, indicating for which
%                   eccentricity to compute the mean Zernike coefficient
%
% Outputs:
%    zCoefs      - Mean Zernike coefficient for a requested J value, eye and
%                   eccentricity
%
% Optional key/value pairs:
%    None.

% Example:
% [wvf, oi] = wvfLoadJaekenArtal2012Data('jIndex', 0:14, 'whichEye','left', 'eccentricity',4)

% See also coneDensityReadData

% First version of function: 04/06/2017 by EK (NYU)

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species','human', @ischar);
p.addParameter('wvfZcoefsSource','JaekenArtal2012',@(x) (ischar(x)));
p.addParameter('jIndex', 0:14, @isnumeric);
p.addParameter('whichEye','left',@(x)(ismember(x,{'left','right'})));
p.addParameter('eccentricity',0, @isnumeric);
p.addParameter('eccentricityUnits','deg',@ischar);
p.addParameter('polarangle',0, @isnumeric);
p.addParameter('polarangleUnits','rad',@ischar);
p.parse(varargin{:});

% Set up params return.
params = p.Results;


%% Load data

% Set species
switch (params.species)
    case {'human'}
        % Set zernike coefficient source
        switch (params.wvfZcoefsSource)
            case {'JaekenArtal2012'}
                % Load the digitized human wvf zernike coefficients from the ISETBio style mat file.  The
                % data file is one matrix, where the first column is
                % patient number, the second column defines zernike coefficient,
                % and the following columns are the measured
                % These each have fields 'density' as a function of 'eccMM' in units of cones/mm2.
                data = rawDataReadData('zCoefsJaekenArtal2012','datatype','isetbiomatfileonpath');
                data = data.data;
            otherwise
                error('Unsupported source specified');
        end
        
    otherwise
        error('Unsupported species specified');
        
end


%% Set parameters to reshape dataset
totalZCoefs         = length(0:14);
totalPatients       = 130;
totalEyes           = length({'left','right'});
totalEccen          = length(-40:1:40);

% Truncate headers and reshape data
data = data(2:end,4:end);
data = reshape(data, totalZCoefs, totalPatients, totalEyes, totalEccen); % zernike x subject x eye x eccentricity

% Preallocate matrix for data
theseZCoef   = wvfOSAIndexToVectorIndex(params.jIndex);
eyeIdx       = strcmp(params.whichEye, {'left','right'});

switch params.polarangle
    case 0
        eccenIdx = find(ismember(-40:1:40, round(params.eccentricity)));
         
    case pi
        eccenIdx = find(ismember(-40:1:40, round(params.eccentricity)));
    otherwise
        error('(wvfLoadJaekenArtal2012Data): Zernike coefficients only available for horizontal meridian')
        return
end

% Truncate data more, to speed up process
data      = data(theseZCoef, :, eyeIdx, eccenIdx);

% Preallocate matrices for psf and otf. We set psf with zernike coefficients
% then convert to OTF and average across patients in OTF space
otfSupportLength = 201;

psf = zeros(otfSupportLength, otfSupportLength, totalPatients);
otf = zeros(otfSupportLength, otfSupportLength, totalPatients);

%% For each requested eccentricity (in degrees) (max 81)
for p = 1:totalPatients
    
    % Create human default wvf
    wvf = wvfCreate;
    
    % Set pupil size to the one reported in Jaeken and Artal (2012)
    wvf = wvfSet(wvf, 'measuredpupilsize', 4);
    wvf = wvfSet(wvf, 'calcpupilsize', 4);
    
    % Set all 15 zernike coeffs at once
    nanidx = isnan(data(:,p));
    if any(nanidx)
        data(nanidx,p) = 0;
    end
    
    wvf =  wvfSet(wvf,'zcoeffs',data(:,p), params.jIndex);
    
    % compute the PSF
    wvf = wvfComputePSF(wvf);
    psf(:,:,p) = wvfGet(wvf, 'psf');
    
    % Get psf support
    psfSupport = wvfGet(wvf, 'spatial Support', 'min');
    
    % % Plot observers PSF
%     figure; subplot(121)
%     xl = [min(psfSupport) max(psfSupport)];  
%     imagesc(xl,xl, psf(:,:,p)); axis square; colormap gray
%     xlabel('X (minutes)')
%     ylabel('Y (minutes)')    
%     title(sprintf('PSF for eccen: %2.0f deg, polar angle: %2.0f rad', round(params.eccentricity), params.polarangle))
    
    % Get otf
    otfSupport = wvfGet(wvf, 'otfSupport', 'deg');
    
    % Get the the OTF
    otf(:,:,p) = wvfGet(wvf,'otf');
    
    % % Plot observers OTF
%     subplot(122)
%     xl = [min(otfSupport) max(otfSupport)];
%     imagesc(xl,xl, abs(otf(:,:,p))); axis square; colormap gray
%     xlabel('X (c/d)')
%     ylabel('Y (c/d)')
%     title(sprintf('OTF for eccen: %2.0f deg, polar angle: %2.0f rad', round(params.eccentricity), params.polarangle))
end


% Take mean of OTF across observers
otfMean = mean(otf, 3);

% Take the amplitudes, shift the frequencies to be centered around zero,
otfMeanAbs = abs(fftshift(otfMean));

% Set the wvf fields to correct OTF and PSF
wvf.otf = {otfMeanAbs};
wvf.psf = {wvfGet(wvf, 'psf')};

% Remove zcoeffs, since the last zcoeffs are from one observer's eye, for
% one eccentricity
wvf.zcoeffs = [];

% Change name
wvf.name = 'Average JaekenArtal WVF';

% Create OI from our wvf
oi = wvf2oi(wvf);

% Get optics
optics = oiGet(oi,'optics');

% Not sure why, but we need to reset the otf data in oi.optics
oi.optics = opticsSet(optics,'otf data',otfMean);
oi.optics.model = 'custom';

% Plot mean of OTF for requested eye
otfSupport = wvfGet(wvf, 'otfSupport', 'deg');
figure; 
subplot(121); 
mesh(otfSupport, otfSupport, wvf.otf{:}); title('OTF of average data Jaeken&Artal2012')

subplot(122);
psfSupport = wvfGet(wvf, 'spatial support', 'min'); 
mesh(psfSupport, psfSupport, wvf.psf{:}); title('PSF of average data Jaeken&Artal2012')






