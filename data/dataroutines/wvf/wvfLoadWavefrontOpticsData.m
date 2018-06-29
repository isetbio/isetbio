function [wvf, oi] = wvfLoadWavefrontOpticsData(varargin)
% Load wave front data from literature, so far only Jaeken & Artal (2012)
% containing Zernike coefficients along the central 80 degrees along horizontal merdian 
%
% Syntax:
%    [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex',[0:14], 'whichEye', ...
%                'left', 'eccentricity', 0, 'whichGroup', 'emmetropes');
%
% Description:
%    Convert a dataset of individual subject Zernike coefficients
%    to a mean OTF for eccentricities -40 to 40 degrees on the horizontal axis
%    for a given eye. Pupil size was 4 mm. Nasal side, 10-18 degrees
%    eccentricity are removed since they are affected by the optic disk.
%       [NOTE: Authors state "severe distorted Hartmann Shack-spots in outermost
%       eccentricities (>35 degrees) are removed from analysis", but don't
%       give a definition of what is considered severe]
%
%    Data are published in Jaeken, B. & Artal, P. (2012) Optical Quality of
%    Emmetropic and Myopic Eyes in the Periphery Measured with High-Angular
%    Resolution. Investigative Ophthalmology & Visual Science, June 2012,
%    Vol. 53, No. 7. Link: https://www.ncbi.nlm.nih.gov/pubmed/22511633
%
%    For information on Zernike coefficient and their names:
%    http://www.telescope-optics.net/monochromatic_eye_aberrations.htm

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
%      None.
%
% Outputs:
%    wvf          - Wavefront aberrations structure based on mean across
%                    subjects in Jaeken & Artal 2012 data
%    oi           - Optical image of off mean across
%                    subjects in Jaeken & Artal 2012 data
%
% Optional key/value pairs:
%    'species'     - Species to estimate for (default 'human')
%                    Options:
%                      'human' - Human.
%    'wvfZcoefsSource' - Data source for coefficients (default
%                    'JaekenArtal2012').
%                    Options:
%                      'JaekenArtal2012' - Jaeken and Artal, 2012. 
%    'jIndex'      - List of OSA J values (default 0:14)
%    'whichEye'    - String defining whether to use right (1) or left (2) eye
%                    (default 'right').
%    'eccentricity' - Integer between -40 and 40 (deg), indicating for which
%                    eccentricity to compute the mean Zernike coefficients for.
%                    (default, 0).
%                      For right eye -40 = nasal, 40 = temporal
%                      For left eye  -40 = temporal, 40 = nasal
%    'eccentricityUnits' - Units in which eccentricity is specified
%                    (default 'deg')
%                    Options:
%                      'deg' - Degrees of visual angle.
%    'whichGroup'  - String or scalar defining which subset of subjects to analyze
%                    (default 'emmetropes').
%                    Options:
%                      'emmetropes' - Emmetropes
%                      'myopes'     - Myopes
%                      'singleRandomEmmetrope' - select a random observer
%                      from the emmetrope group
%                      single scalar between 1-130 - select a specific observer from the entire pool 
%    'verbose'     - Boolean (default false). Print things out.
%
% Examples are included in the source code.
%
% See also coneDensityReadData
%

% Examples:
%{
    [wvf, oi] = wvfLoadWavefrontOpticsData('source', 'JaekenArtal2012', 'jIndex', 0:14, ...
                'whichEye', 'left', 'eccentricity', 4, 'whichGroup', 'emmetropes', ...
                'verbose', true);
%}

% History:
%   04/06/18    ek (NYU) First version of function
%   05/05/18    dhb      Cosmetic.
%
%% 0. Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species','human',@ischar);
p.addParameter('wvfZcoefsSource','JaekenArtal2012',@ischar);
p.addParameter('jIndex', 0:14, @isnumeric);
p.addParameter('whichEye','right',@(x)(ismember(x,{'right','left'})));
p.addParameter('eccentricity',0, @isnumeric);
p.addParameter('eccentricityUnits','deg',@ischar);
p.addParameter('whichGroup','emmetropes',@(x) (isscalar(x) | ischar(x)) );
p.addParameter('verbose',false,@islogical);
p.parse(varargin{:});

% Set up params return.
params = p.Results;

%% 1. Load data
%
% Set units
switch (params.eccentricityUnits)
    case 'deg'
    otherwise
        error('Unsupported units specified');
end

% Set species
switch (params.species)
    case {'human'}
        % Set zernike coefficient source
        switch (params.wvfZcoefsSource)
            case {'JaekenArtal2012'}
                % Load the digitized human wvf zernike coefficients from the ISETBio style mat file.  The
                % data file is one matrix, where the first column is
                % subject number, the second column defines zernike coefficient,
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

%% 2. Set parameters to reshape dataset
totalZCoefs         = length(0:14);
totalSubjects       = 130;
totalEyes           = length({'right','left'});
totalEccen          = length(-40:1:40);

% Truncate headers and reshape data
data = data(2:end,4:end);
data = reshape(data, totalZCoefs, totalSubjects, totalEyes, totalEccen); % zernike x subject x eye x eccentricity

% Preallocate matrix for data
theseZCoef   = wvfOSAIndexToVectorIndex(params.jIndex);
eyeIdx       = strcmp(params.whichEye, {'right','left'});
eccenIdx     = ismember(-40:1:40, round(params.eccentricity));

% Remove optic disk for right (1) and left (2) eye
opticDisk.RE = ismember(-40:40, -18:-10);
opticDisk.LE = ismember(-40:40, 10:18);
data(:,:,1, opticDisk.RE) = NaN;
data(:,:,2, opticDisk.LE) = NaN;

% Report requested eccentricities and retinal side
switch params.whichEye
    case 'right' % -40:40 corresponds to nasal to temporal retina
        % Check eccentricities for optic disk
        if ismember(-18:-10, round(params.eccentricity)); warning('Some requested eccentricities fall within optic disk\n'); end
        nas = round(params.eccentricity) <=0;
        temp  = round(params.eccentricity) >=0;
        if (params.verbose)
            fprintf('\nData are given for right eye, nasal retina, eccen: %d \t temporal retina, eccen: %d\n', params.eccentricity(nas), params.eccentricity(temp));
        end
            
    case 'left' % -40:40 corresponds to temporal to nasal retina
        % Check eccentricities for optic disk
        if ismember(10:18, round(params.eccentricity)); warning('Some requested eccentricities fall within optic disk'); end
        temp = round(params.eccentricity) <=0;
        nas  = round(params.eccentricity) >=0;
        if (params.verbose)
            fprintf('\nData are given for left eye, temporal retina, eccen: %d \t nasal retina, eccen: %d\n', params.eccentricity(temp), params.eccentricity(nas));
        end
        
    otherwise
        error('Unknown eye specified.  Use ''right'' or ''left''');
end

% Select subject group
[~, emmetropes, myopes, group] = wvfSortSubjectDataJaekenArtal2012('verbose',params.verbose);
switch params.whichGroup
    case 'emmetropes'
        if find(eyeIdx) == 1; subjectIdx = getfield(emmetropes, 'RE');
        else subjectIdx = getfield(emmetropes, 'LE'); end
    case 'myopes'
        if find(eyeIdx) == 1; subjectIdx = getfield(myopes, 'RE');
        else subjectIdx = getfield(myopes, 'LE'); end
    case 'singleRandomEmmetrope'
        if find(eyeIdx) == 1; subjects = getfield(emmetropes, 'RE');
        else subjects = getfield(emmetropes, 'LE'); end        
        selectOne = randi(length(subjectIdx));
        subjectIdx = subjectIdx(selectOne(1));        
    otherwise
        if isscalar(params.whichGroup)
            subjectIdx = params.whichGroup;
            fprintf('\n Subject index: %d\n', subjectIdx)
        else    
            error('Unknown group specified');
        end
end

% Truncate data to speed up process
data      = data(theseZCoef, subjectIdx, eyeIdx, eccenIdx);

% Preallocate matrices for psf and otf. We set psf with zernike coefficients
% then convert to OTF and average across subjects in OTF space
otfSupportLength = 201;
psf = zeros(otfSupportLength, otfSupportLength, totalSubjects);
otf = zeros(otfSupportLength, otfSupportLength, totalSubjects);

%% 3. For each requested eccentricity (in degrees) (max 81)
for p = 1:length(subjectIdx)
    
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
    wvf =  wvfSet(wvf, 'measured wave', 780);

    % compute the PSF
    wvf = wvfComputePSF(wvf);
    psf(:,:,p) = wvfGet(wvf, 'psf');

    % Optional plots for debugging
    % Get psf support
    % psfSupport = wvfGet(wvf, 'spatial Support', 'deg');

    % % Plot observers PSF

    % figure; subplot(121)
    % xl = [min(psfSupport) max(psfSupport)];  
    % imagesc(xl,xl, psf(:,:,p)); axis square; colormap gray
    % xlabel('X (minutes)')
    % ylabel('Y (minutes)')    
    % title(sprintf('PSF for eccen: %2.0f deg, polar angle: %2.0f rad', round(params.eccentricity), params.polarangle))   

    % Get otf support
    % otfSupport = wvfGet(wvf, 'otfSupport', 'deg');
    
    % Get the  OTF
    otf(:,:,p) = wvfGet(wvf,'otf');

    % % Plot observers OTF
    % subplot(122)
    % xl = [min(otfSupport) max(otfSupport)];
    % imagesc(xl,xl, abs(otf(:,:,p))); axis square; colormap gray
    % xlabel('X (c/d)')
    % ylabel('Y (c/d)')
    % title(sprintf('OTF for eccen: %2.0f deg, polar angle: %2.0f rad', round(params.eccentricity), params.polarangle))

end 

% Take mean of OTF across observers
otfMean = nanmean(otf, 3);

% Take the amplitudes, (used to also use fftshit, to shift the frequencies
% to be centered around zero)
otfMeanAbs = abs((otfMean));

% Set the wvf fields to correct OTF and PSF
wvf.otf = {otfMeanAbs};
wvf.psf = {wvfGet(wvf, 'psf')};

% Remove zcoeffs, since the last zcoeffs are from one observer's eye, for
% one eccentricity
wvf.zcoeffs = [];

% Change name
wvf.name = 'Average Subject JaekenArtal WVF';

% Set wavelength
wvf = wvfSet(wvf,'calc wave',550);

% Create OI from our wvf
oi = wvf2oi(wvf);

% Get optics
optics = oiGet(oi,'optics');

% [EK:] Not sure why, but we need to reset the otf data in oi.optics
oi.optics = opticsSet(optics,'otf data',otfMeanAbs);
oi.optics.model = 'custom';

% Plot OTF and PSF
%{
if (params.verbose)
    wvfPlot(wvf,'2d otf','um',550); xlim([-0.6 0.6]), ylim([-0.6 0.6])

    % Plot PSF based off the OTF mean requested eye
    wvfPlot(wvf,'2d psf space','deg',550);
end
%}

return


