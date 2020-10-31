function [wvf, oi] = wvfLoadWavefrontOpticsData(varargin)
% Load wavefront data from literature for Zernike coeffs
%
% Syntax:
%   [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', [0:14], ...
%       'whichEye', 'left', 'eccentricity', 0, 'whichGroup', 'emmetropes');
%
% Description:
%    This function loads the wavefront data from literature, including
%    Jaeken & Artal (2012) and Polans et al. (2015) containing Zernike
%    coefficients along the central 80 degrees along horizontal meridian
%    and central 50 degrees along the vertical meridian.
%
%    Convert a dataset of individual subject Zernike coefficients
%    to a mean OTF for the requested eccentriciy for a given eye.
%
%    Jaeken dataset:
%       - Pupil size was 4 mm.
%       - Zernike coefficients Z0-Z14
%       - Data from -40:1:40 degrees on the horizontal axis
%       - Contains data from left and right eye
%       - Right eye, positive / negative = temporal / nasal
%         Left eye, positive/ negative = nasal / temporal
%       - Nasal side, 10-18 degrees eccentricity will be removed since they
%         are affected by the optic disk.
%
%       [Note: Authors state "severe distorted Hartmann Shack-spots in
%       outermost eccentricities (>35 degrees) are removed from analysis",
%       but don't give a definition of what is considered severe]
%
%    Polans dataset:
%       - Pupil size was 4 mm.
%       - Zernike coefficients Z3-Z20.
%       - Z4 has -0.8D of chromatic aberations taken into account
%       - Data from -40:1:40 degrees on the horizontal axis and -25:5:25
%         on the vertical axis.
%       - Horizontal negative: nasal retina
%         Horizontal positive: temporal retina

%         Vertical negative: inferior retina
%         Vertical positive: superior retina
%       - Contains data from right eye only.
%       - Nasal side, 10-18 degrees eccentricity will be removed since they
%         are affected by the optic disk.
%
%       [Note: Authors state "severe distortions outermost
%       eccentricities (>35 degrees) are removed from analysis", but don't
%       give a definition of what is considered severe]
%
%    Please note that there are examples included in the source code below.
%    To access these examples, type 'edit wvfLoadWavefrontOpticsData.m'
%    into the Command Window.
%
% Inputs:
%    None.
%
% Outputs:
%    wvf                  - Struct. Wavefront aberrations structure based
%                           on the mean across subjects in Jaeken & Artal
%                           2012 data.
%    oi                   - Struct. Optical image of off mean across
%                           subjects in Jaeken & Artal 2012 data.
%
% Optional key/value pairs:
%    'species'            - String. Species to estimate for. (Default
%                           'human') Options include:
%               'human': Human.
%    'wvfZcoefsSource'    - String. Data source for coefficients (default
%                           'JaekenArtal2012'). Options include:
%               'JaekenArtal2012': Jaeken and Artal, 2012.
%               'Polans2015': Polans et al., 2015.
%    'jIndex'             - Vector. List of OSA J values (default 0:14)
%    'whichEye'           - String. String defining whether to use 'right'
%                           or 'left' eye (default 'right').
%    'eccentricity'       - Vector. A 1x2 vector of integers. The first
%                           integer represents the horizontal eccentricity,
%                           and the second the vertical eccentricity.  If
%                           only a single value is passed, it is assumed to
%                           be horizontal eccentricity and vertical
%                           eccentricity is taken to be 0.
%                           (Default is [0, 0]).
%               HORIZONTAL: Integer between -40 and 40 (deg), steps of 1, 
%                           indicating for which eccentricity to compute
%                           the mean Zernike coefficients for on the
%                           horizontal axis. For the right eye -40 =
%                           temporal, 40 = nasal. For the left eye  -40 =
%                           nasal, 40 = temporal.
%               VERTICAL:   Integer between -25 and 25 (deg), steps of 5,
%                           indicating for which eccentricity to compute
%                           the mean Zernike coefficients for on the
%                           vertical axis. For the right eye only: -25 =
%                           inferior, 25 = superior.
%    'eccentricityUnits'  - Units in which eccentricity is specified
%                           (default 'deg') Options include:
%               'deg':        Degrees of visual angle.
%    'whichGroup'         - String or scalar defining which subset of the
%                           subjects to analyze (default 'emmetropes'). The
%                           options include:
%               'emmetropes': Emmetropes
%               'myopes': Myopes
%               'singleRandomEmmetrope': Select a random observer from the
%                   emmetrope group single scalar between 1-130 (Jaeken)
%               <scalar value>: Select a specific observer from
%                   the entire pool (1-130, Jaeken) or 1-10 (Polans).
%    'relativeRefraction' - Boolean. A boolean to state whether you want to
%                           compute the relative refraction from the fovea.
%                           If true, then the foveal Z4 (Defocus) component
%                           will be subtracted from all Z4 components.
%    'verbose'            - Boolean. Default false. Print things out.
%
% References:
%    Jaeken, B. & Artal, P. (2012) Optical Quality of
%    Emmetropic and Myopic Eyes in the Periphery Measured with High-Angular
%    Resolution. Investigative Ophthalmology & Visual Science, June 2012, 
%    Vol. 53, No. 7. Link: https://www.ncbi.nlm.nih.gov/pubmed/22511633
%
%    Polans J, Jaeken, B., McNabb, R.P., Artal, P., Izatt, J.A. (2015)
%    Wide-field optical model of the human eye with asymmetrically tilted
%    and decentered lens that reproduces measured ocular aberrations.
%    Optica, Vol. 2, No. 2. Link: https://doi.org/10.1364/OPTICA.2.000124
%
%    For information on Zernike coefficient and their names:
%    http://www.telescope-optics.net/monochromatic_eye_aberrations.htm
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
% See Also:
%    coneDensityReadData
%

% History:
%    04/06/18  ek (NYU) First version of function
%    05/05/18  dhb      Cosmetic.
%    09/26/18  jnm      Formatting. Add catch for single eccentricity
%                       value, note that examples are BROKEN, added TODO.
%    10/19/20  dhb      Added comment that if only a single eccentricity is
%                       passed, it is taken to be the horizontal eccentricity with vertical
%                       eccentricity as 0.  Removed warning for this case,
%                       but print a message if verbose is true.
%                       Fixed examples not to pass a single number for eccentricity.

% Examples:
%{
    [wvf, oi] = wvfLoadWavefrontOpticsData(...
        'wvfZcoefsSource', 'JaekenArtal2012', 'jIndex', 0:14, ...
        'whichEye', 'left', 'eccentricity', [4 0], ...
        'whichGroup', 'emmetropes', 'verbose', true);
%}
%{
    [wvf, oi] = wvfLoadWavefrontOpticsData(...
        'wvfZcoefsSource', 'Polans2015', 'jIndex', 3:14, ...
        'whichEye', 'right', 'eccentricity', [5, 5], ...
        'whichGroup', 1:10, 'relativeRefraction', true, 'verbose', true);
%}



%% 0. Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species', 'human', @ischar);
p.addParameter('wvfZcoefsSource', 'JaekenArtal2012', @ischar);
p.addParameter('jIndex', 0:14, @isnumeric);
p.addParameter('whichEye', 'right', @(x)(ismember(x, {'right', 'left'})));
p.addParameter('eccentricity', [0, 0], @isnumeric);
p.addParameter('eccentricityUnits', 'deg', @ischar);
p.addParameter('whichGroup', 'emmetropes', ...
    @(x) (isscalar(x) | ischar(x)) | ismatrix(x));
p.addParameter('relativeRefraction', false, @islogical);
p.addParameter('verbose', false, @islogical);
p.parse(varargin{:});

% Set up params return.
params = p.Results;

%% 1. Load data
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
                % Load the digitized human wvf zernike coefficients from
                % the ISETBio style mat file.  The data file is one matrix,
                % where the first column is subject number, the second
                % column defines zernike coefficient, and the following
                % columns are the measured
                allData = rawDataReadData('zCoefsJaekenArtal2012', ...
                    'datatype', 'isetbiomatfileonpath');
                allData = allData.data;

                horzEccen = allData(1, 4:end);
                vertEccen = 0;
                usedZCoefs = 0:14;
                totalZCoefs = length(usedZCoefs);
                totalSubjects = sum(~isnan(unique(allData(:, 1))));
                totalEyes = length({'right', 'left'});
                totalEccen = length(horzEccen);

                % Truncate headers and reshape data
                allData = allData(2:end, 4:end);
                % zernike x subject x eye x eccentricity
                allData = reshape(allData, totalZCoefs, totalSubjects, ...
                    totalEyes, totalEccen);

                % Get central refractive error of each subject if requested
                if params.relativeRefraction
                    eyeIdx = strcmp(params.whichEye, {'right', 'left'});
                    Z4idx = 5; % since zcoeffs are from Z0-Z14
                    centralRefraction = ...
                        allData(Z4idx, :, eyeIdx, horzEccen == 0);
                end

                % Remove optic disk for right (1) and left (2) eye
                opticDisk.RE = ismember(horzEccen, -18:-10);
                opticDisk.LE = ismember(horzEccen, 10:18);
                allData(:, :, 1, opticDisk.RE) = NaN;
                allData(:, :, 2, opticDisk.LE) = NaN;

            case {'Polans2015'}
                % Load the digitized human wvf zernike coefficients from
                % the ISETBio style mat file.  The data file is one matrix,
                % where the first column is subject number, the second
                % column defines measured eccentricity (-25:5:25 vertical x
                % -40:1:40 horizontal) and the following columns are the
                % measured zernike coefficients from ANSI2005 Z4:Z20

                allData = rawDataReadData('zCoefsPolans2015', ...
                    'datatype', 'isetbiomatfileonpath');
                allData = allData.data;

                % Note, flipped horizontal axis from Jaeken and Artal data
                horzEccen = allData(1, 1:81, 2);
                vertEccen = allData(1, 1:81:end, 1);
                usedZCoefs = 3:20;
                totalZCoefs = length(usedZCoefs);
                totalSubjects = size(allData, 1);
                totalEyes = length({'right'});
                totalEccenHorz = length(horzEccen);
                totalEccenVert = length(vertEccen);
                
                % Truncate coordinate headers
                allData = allData(:, :, 3:end);
                
                % Permute so that data matrix goes from
                %   OLD (subject x all coords x zernike) to
                %   NEW (zernike, subject, horz coords x vert coords)

                % zernike, subject, vertical*horizontal
                allData = permute(allData, [3, 1, 2]);
                % zernike x subject x vertical x horizontal
                allData = reshape(allData, totalZCoefs, totalSubjects, ...
                    totalEccenHorz, totalEccenVert);
                % zernike x subject x horizontal x vertical
                % allData = permute(allData, [1, 2, 4, 3]);

                % Check requested eye, can only be right
                assert(strcmp(params.whichEye, 'right'));
                
                % Get central refractive error of each subject if requested
                if params.relativeRefraction
                    Z4idx = 2; % since zcoeffs are from Z3-Z20
                    centralRefraction = allData(Z4idx, :, ...
                        horzEccen == 0, vertEccen == 0);
                end

                % Remove optic disk for right (1), since only data from
                % right eye
                opticDisk.RE = ismember(horzEccen, -18:-10);
                allData(:, :, opticDisk.RE, :) = NaN;
                
            otherwise
                error('Unsupported source specified');
        end

    otherwise
        error('Unsupported species specified');
end

%% 2. Use parameters to reshape dataset
% Preallocate matrix for data
theseZCoef = wvfOSAIndexToVectorIndex(params.jIndex);
eyeIdx = strcmp(params.whichEye, {'right', 'left'});

% Check if eccentricity is a single value, fill in vertical as zero if so.
if length(params.eccentricity) < 2
    params.eccentricity = [params.eccentricity, 0];
    if (params.verbose)
        fprintf('Only a single value for eccentricity supplied. Assuming vertical eccentricity is 0');
    end
end

% Horizontal followed by Vertical
eccenIdx = [find(ismember(horzEccen, round(params.eccentricity(1)))), ...
                find(ismember(vertEccen, round(params.eccentricity(2))))];
    
%% 3. Report requested eccentricities and retinal side
switch params.whichEye
    case 'right' % -40:40 corresponds to nasal to temporal retina
        % Check eccentricities for optic disk
        if ismember(-18:-10, round(params.eccentricity(1)))
            warning(strcat("Some of the requested eccentricities ", ...
                "fall within optic disk\n"));
        end

        % Check horizontal axis (nasal or temporal)
        nas = round(params.eccentricity(1)) < 0;
        temp = round(params.eccentricity(1)) > 0;

        % Check vertical axis (inferior or superior)
        inf = round(params.eccentricity(2)) < 0;
        sup = round(params.eccentricity(2)) > 0;

        % Report coordinates for requested data:
        if (params.verbose)
            fprintf('\nData are given for right eye:\n')
            if nas
                fprintf('nasal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            elseif temp
                fprintf('temporal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            else
                fprintf('foveal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            end
            
            if inf
                fprintf('inferior retina, eccen: %2.0f\n', ...
                    params.eccentricity(2))
            elseif sup
                fprintf('superior retina, eccen: %2.0f\n', ...
                    params.eccentricity(2))
            else
                fprintf('foveal retina, eccen: %2.0f\n', ...
                    params.eccentricity(2))
            end
        end
        
    case 'left' % -40:40 corresponds to temporal to nasal retina
        % Check eccentricities for optic disk
        if ismember(10:18, round(params.eccentricity(1)))
            warning(strcat("Some requested eccentricities fall ", ...
                "within optic disk"));
        end

        % Check horizontal axis (nasal or temporal)
        temp = round(params.eccentricity(1)) <0;
        nas = round(params.eccentricity(1)) <0;

        % Report coordinates for requested data:
        if (params.verbose)
            fprintf('\nData are given for left eye:\n')
            if nas
                fprintf('nasal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            elseif temp
                fprintf('temporal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            else
                fprintf('foveal retina, eccen: %2.0f\n', ...
                    params.eccentricity(1))
            end
        end

    otherwise
        error("Unknown eye specified.  Use 'right' or 'left'");
end

%% 4. Report requested subject group or single subject number from data
if strcmp(params.wvfZcoefsSource, 'JaekenArtal2012')
    % Select subject group
    [~, emmetropes, myopes, ~] = ...
        wvfSortSubjectDataJaekenArtal2012('verbose', params.verbose);
    switch params.whichGroup
        case 'emmetropes'
            if find(eyeIdx) == 1, subjectIdx = emmetropes.RE;
            else, subjectIdx = emmetropes.LE; end
        case 'myopes'
            if find(eyeIdx) == 1, subjectIdx = myopes.RE;
            else, subjectIdx = myopes.LE; end
        case 'singleRandomEmmetrope'
            if find(eyeIdx) == 1, subjects = emmetropes.RE;
            else, subjects = emmetropes.LE; end
            selectOne = randi(length(subjects));
            subjectIdx = subjects(selectOne(1));
        otherwise % single subject
            if isscalar(params.whichGroup)
                subjectIdx = params.whichGroup;
                fprintf('\n Subject index: %d\n', subjectIdx)
            else
                error('Unknown group specified');
            end
    end

    % Truncate data to speed up process
    currData = allData(theseZCoef, subjectIdx, eyeIdx, eccenIdx(1));
    if isnan(currData)
        fprintf('Clean up the data. ecc %d\n',params.eccentricity(1));
        error('NaNs for subject.'); 
    end

elseif strcmp(params.wvfZcoefsSource, 'Polans2015')
    % Only option is to get single subject
    if isscalar(params.whichGroup) || ismatrix(params.whichGroup)
        subjectIdx = params.whichGroup;
        fprintf('\n Subject index: %d\n', subjectIdx)
    else
        error('Unknown subject specified');
    end

    % Truncate data to speed up process
    theseZCoef = find(ismember(params.jIndex, usedZCoefs));
    currData = allData(theseZCoef, subjectIdx, eccenIdx(1), eccenIdx(2));
end

%% 5. Preallocate matrices for psf and otf. 
otfSupportLength = 201;
all_psf = zeros(otfSupportLength, otfSupportLength, length(subjectIdx));
all_otf = zeros(otfSupportLength, otfSupportLength, length(subjectIdx));
usedSubjectData = zeros(length(theseZCoef), length(subjectIdx));
subjectWithoutData = zeros(length(subjectIdx), 1);

%% 3. Get psf with zernike coefficients then convert to OTF 
% for the requested coordinates (in degrees)
for p = 1:length(subjectIdx)
    % Get subject data
    subjectZData = currData(:, p);

    % If requested to use a relative measure of wavefront, i.e. use the
    % foveal refraction as the baseline, then subtract the refractive error
    % of foveal Z4 from all other Z4 eccentricities.
    if params.relativeRefraction
        subjectCRF = centralRefraction(p);
        fprintf('Subject foveal refraction: %1.2f\n', subjectCRF)
        if any(ismember(params.jIndex, 4))
            idx = find(ismember(params.jIndex, 4));
            subjectZData(idx) = subjectZData(idx) - subjectCRF;
        end
    end

    % Check for all nans or all zero's, if so, do not use these data
    nanidx = any(isnan(subjectZData));
    zeroidx = (sum(subjectZData) == 0);
    if any(nanidx) || any(zeroidx)
        subjectWithoutData(p) = 1;
    else
        % Create human diffraction limited PSF for 550nm light & 3mm pupil
        wvf = wvfCreate;
        
        % Set pupil size to the one reported in both Jaeken and Artal
        % (2012) and Polans (2015) (i.e. 4 mm) and wavelength 550nm.
        wvf = wvfSet(wvf, 'measuredpupilsize', 4);
        wvf = wvfSet(wvf, 'calcpupilsize', 4);
        wvf = wvfSet(wvf, 'calc wave', 550);

        % Set all zernike coeffs at once
        wvf = wvfSet(wvf, 'zcoeffs', subjectZData, params.jIndex);

        % Report the zernikes if requested
        if params.verbose
            if length(subjectIdx)<2
            fprintf('Zernike coefficients (um):\n')
            disp(subjectZData)
            end
        end

        % Compute and accumulate the PSF, OTF and zernikes across subjects
        wvf = wvfComputePSF(wvf);
        all_psf(:, :, p) = wvfGet(wvf, 'psf');
        all_otf(:, :, p) = wvfGet(wvf, 'otf');
        usedSubjectData(:, p) = subjectZData;

        % Optional plots for debugging
%{
if (params.verbose)
        % Plot observers PSF

        psfSupport = wvfGet(wvf, 'psfangularsamples', 'deg');
        figure(1);
        clf;
        subplot(121)
        xl = [min(psfSupport) max(psfSupport)];
        imagesc(xl, xl, all_psf(:, :, p == subjectIdx));
        axis square;
        colormap gray
        xlabel('X (minutes)')
        ylabel('Y (minutes)')
        title(sprintf('PSF for horz eccen: %2.0f deg, vert: %2.0f deg', ...
            round(params.eccentricity(1)), params.eccentricity(2)))

        % Get otf support
        otfSupport = wvfGet(wvf, 'otfSupport', 'deg');

        % Plot observers OTF
        subplot(122)
        xl = [min(otfSupport) max(otfSupport)];
        imagesc(xl, xl, abs(all_otf(:, :, p)));
        axis square;
        colormap gray
        xlabel('X (c/d)')
        ylabel('Y (c/d)')
        title(sprintf('OTF for horz eccen: %2.0f deg, vert: %2.0f deg', ...
            round(params.eccentricity(1)), params.eccentricity(2)))
        end
%}        
    end
    
end

%% If more than one subject requested, average across subjects in OTF space
if length(subjectIdx) > 1
    % First check for subjects without data
    if any(subjectWithoutData)
        fprintf('(%s): Excluding subject nr %02d because there is no data available \n', ...
            mfilename, find(subjectWithoutData))
    end
    all_otf = all_otf(:, :, ~subjectWithoutData);

    % Take mean of OTF across observers
    otfMean = nanmean(all_otf, 3);

    % Take the amplitudes, (used to also use fftshit, to shift the
    % frequencies to be centered around zero)
    otfMeanAbs = abs(otfMean);

    % Set the wvf fields to correct OTF and PSF
    wvf.otf = {otfMeanAbs};
    wvf = wvfComputePSF(wvf);
    wvf.psf = {wvfGet(wvf, 'psf')};

    % Change name
    wvf.name = sprintf('Average Subject Imported WVF, %s', ...
        params.wvfZcoefsSource);
    
    % Add mean zcoeffs to wavefront structure
    wvf.zcoeffs = mean(usedSubjectData(:, ~subjectWithoutData), 2);
else % No averaging, just use single subject, add otf and change name

    % Add otf
    wvf.otf = {wvfGet(wvf, 'otf')};

    % Change name
    wvf.name = sprintf('Single Subject Imported WVF, %s', ...
        params.wvfZcoefsSource);
end

% Create OI from our wvf
oi = wvf2oi(wvf);

% Get optics
optics = oiGet(oi, 'optics');

% [Note: EK - Not sure why, but we need to reset the otf data in oi.optics]
oi.optics = opticsSet(optics, 'otf data', wvf.otf{1});
oi.optics.model = 'custom';

% Plot OTF and PSF
%{
if (params.verbose)
    wvfPlot(wvf, '2d otf', 'um', 550);
    xlim([-0.6 0.6]);
    ylim([-0.6 0.6]);

    % Plot PSF based off the OTF mean requested eye
    wvfPlot(wvf, '2d psf angle', 'deg', 550);
end
%}

return
