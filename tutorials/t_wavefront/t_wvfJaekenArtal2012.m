%% t_wvfJaekenArtal2012
%
% Description:
%    Tutorial showing average point spread function (PSF) and average optical
%    transfer function (OTF) for emmetropes and myopes based on Jaeken & Artal 
%    2012 dataset.
%
%    In short, Jaeken and Artal dataset provides higher order aberrations
%    along the horizontal meridian (central 80 degrees, sampled at 1 degree)
%    for both eyes. These data contain 15 zernike coefficients (OSA j-indices:
%    0:14) for each sample, for 130 subjects. Data are based on a 4 mm pupil
%    diameter, measured at 780 nm. No correction for chromatic aberration.
%
%    The subjects can be divided into emmetropes based on their mean central
%    refraction in diopters (central 5 degrees, using the defocus only (OSA
%    j-index = 4). In the corresponding published article, the subjects are
%    also dividied into 6 non-overlapping groups based on the strength of
%    central refractions. This division is visualized by the function
%    sortPatientDataJaekenArtal2012 (inside wvfLoadJaekenArtal2012Data).
%
%    Reference:
%       Jaeken, B. & Artal, P. (2012) Optical Quality of Emmetropic and Myopic
%       Eyes in the Periphery Measured with High-Angular Resolution. Investigative
%       Ophthalmology & Visual Science, June 2012, Vol. 53, No. 7
%       Link: https://www.ncbi.nlm.nih.gov/pubmed/22511633
%
% See also: wvfLoadJaekenArtal2012Data and wvfSortSubjectDataJaekenArtal2012
%

% History:
%   05/03/18    ek (NYU) First version of function
%   05/05/18    dhb      Cosmetic.

%% Define which zernike coefficients we want to use:
zIndices    = 0:14;    % In this case, all of them
whichEye    = 'left';
eccen       = 4;       % degrees
whichGroup  = 'emmetropes';

%% Get wavefront and optics from Artal data with the requested parameters:
% The function wvfLoadJaekenArtal2012Data loads the wavefront zernike
% aberration data, and reconstructs one PSF per subject, then converted to
% an OTF per subject, then we average the subject's OTFs. Lastly, the
% average OTF will get converted back to an average PSF. The individual
% PSFs are constructed under the measured wavelength (780 nm), but then
% calculated and plotted for a more sensible (i.e. in the range of human
% sensitivity) wavelength (550 nm).
[wvf, oi] = wvfLoadWavefrontOpticsData('source', 'JaekenArtal2012', 'jIndex', zIndices, 'whichEye', whichEye, ...
    'eccentricity', eccen, 'whichGroup', whichGroup, 'verbose', true);

%% Visualize

% Get [x,y] support for plotting OTF
otfSupport = wvfGet(wvf, 'otfSupport', 'mm');

% Plot OTF
vcNewGraphWin; 
surf(otfSupport, otfSupport, fftshift(wvf.otf{1}));
set(gca, 'XLim', [-100 100], 'YLim', [-100 100])
xlabel('Freq (lines/mm)'); ylabel('Freq (lines/mm)');
title(sprintf('%s: OTF 550 nm, pupil 4 mm, eccen %d deg, %s eye', whichGroup, eccen, whichEye))

% Get [x,y] support for plotting PSF
psfSupport  = wvfGet(wvf, 'spatial Support', 'um');
centeredPSF = [wvf.psf{1}(101:end,:); wvf.psf{1}(1:100,:)];
centeredPSFNormalized = centeredPSF./sum(centeredPSF);

% Plot PSF
vcNewGraphWin; 
surf(psfSupport, psfSupport, centeredPSF)
set(gca, 'XLim', [-40 40], 'YLim', [-40 40])
xlabel('Pos (um)'); ylabel('Pos (um)');
title(sprintf('%s: PSF 550 nm, pupil 4 mm, eccen %d deg, %s eye', whichGroup, eccen, whichEye))

% Plot normalized PSF
vcNewGraphWin; 
surf(psfSupport, psfSupport, centeredPSFNormalized)
set(gca, 'XLim', [-40 40], 'YLim', [-40 40])
xlabel('Pos (um)'); ylabel('Pos (um)');
title(sprintf('%s: Normalized PSF 550 nm, pupil 4 mm, eccen %d deg, %s eye', whichGroup, eccen, whichEye))
