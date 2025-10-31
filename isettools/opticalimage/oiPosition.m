function [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase,varargin)
% OIPOSITION Return oi, psf, and related parameters for a subject at a spatial position
%
% Synopsis
%   [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase,varargin)
%
% Brief description
%   Wrapper function to compute the Optical Image (OI) and Point
%   Spread Function (PSF) for a selected subject from a specified
%   database (Artal2012 or Polans2015), at a given retinal
%   eccentricity and pupil size.
%
%   N.B. The Artal and Polans analyses all assume 290 micronsPerDegree
%   at the retina.
%
% Inputs
%   zCoeffDatabase  - String, must be 'Polans2015' or 'Artal2012'.
%
% Optional Key/Value Pairs
%   'position'        - [H V] vector of eccentricity in degrees. Default is [3 1].
%   'pupildiameter'   - Scalar, pupil diameter in millimeters. Default is 3 mm.
%   'subjectrank'     - Scalar, rank order of the subject based on foveal PSF quality (Strehl).
%                       1 is the best subject available in the database. Default is 1.
%   'wave'            - Scalar, wavelength in nanometers for which to compute the optics.
%                       The underlying function uses this as the single wavelength to compute. Default is 550 nm.
%   'eyeside'         - String, 'left' or 'right' eye. Default is 'right'.
%   'centerpsf'       - Logical, if true, the PSF is centered around its peak. Default is true.
%
% Returns
%   oi              - The ISETBio Optical Image structure.
%   psf             - The computed point-spread function (PSF).
%   support         - Struct containing spatial and wavelength support information:
%                     .psfX, .psfY: PSF support in microns.
%                     .wave: The wavelength used in the computation (nm).
%                     .position: The [H V] eccentricity position used (degrees).
%   zCoeffs         - The interpolated Zernike coefficients (microns) used for the optics.
%   subjID          - The selected subject ID from the database.
%
% See also
%   ArtalOptics.oiForSubjectAtEccentricity, PolansOptics.oiForSubjectAtEccentricity,
%   PolansOptics.constants, ArtalOptics.constants
%
% Example:
% % Get Polans subject with best foveal quality at 5 deg horizontal, 2 mm pupil
% [oi, psf, support] = oiPosition('Polans2015', 'position', [5 0], 'pupildiameter', 2);
% figure; mesh(support.psfX,support.psfY,psf);
%

% Example:
%{
[oi,psf] = oiPosition('Artal2012');
%}
%{
[oi,psf,support,zcoeffs,subject] = oiPosition('Polans2015');
ieFigure;
mesh(support.psfX,support.psfY,psf);
xlabel('Position (um)');
ylabel('Position (um)');
title(sprintf('Subject %d - Wave %d - Position [%d,%d]',subject,support.wave,support.position));
%}

%% Parse parameters

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('zCoeffDatabase', @(x)(ismember(x, {'Polans2015', 'Artal2012'})));
p.addParameter('position',[3 0],@isvector);
p.addParameter('pupildiameter',3,@isnumeric)   % millimeters
p.addParameter('subjectrank',1,@isnumeric);    % PSF quality (Strehl)
p.addParameter('wave',550,@isscalar);          % nanometers
p.addParameter('eyeside','right',@(x)(ismember(x,{'left','right'})));
p.addParameter('centerpsf',true,@islogical);

p.parse(zCoeffDatabase,varargin{:});

wave     = p.Results.wave;
position = p.Results.position;
whichEye = p.Results.eyeside;
pupilDiamMM  = p.Results.pupildiameter;
centerPSF    = p.Results.centerpsf;
subjectrank  = p.Results.subjectrank;

micronsPerDegree = 290;
wavefrontSpatialSamples = 201;


%% Get the data, make the oi
switch (zCoeffDatabase)
    case 'Polans2015'
        switch whichEye
            case 'left'
                whichEye = PolansOptics.constants.leftEye;
            case 'right'
                whichEye = PolansOptics.constants.rightEye;
        end

        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSubjectsID = PolansOptics.constants.subjectRanking;
        subjID = rankedSubjectsID(min([numel(rankedSubjectsID) subjectrank]));

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(subjID);

        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            PolansOptics.oiForSubjectAtEccentricity(subjID, ...
                whichEye, position, pupilDiamMM, wave, micronsPerDegree, ...
               'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
               'subtractCentralRefraction', subtractCentralRefraction, ...
               'zeroCenterPSF', centerPSF);       

    case 'Artal2012'
        switch whichEye
            case 'left'
                whichEye = ArtalOptics.constants.leftEye;
            case 'right'
                whichEye = ArtalOptics.constants.rightEye;
        end

        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSubjectsID = ArtalOptics.constants.subjectRanking(whichEye);
        subjID = rankedSubjectsID(min([numel(rankedSubjectsID) subjectrank]));

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, subjID);
        if (position(2) ~= 0)
            fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
            position(2) = 0;
        end
        
        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            ArtalOptics.oiForSubjectAtEccentricity(subjID, ...
            whichEye, position, pupilDiamMM, wave, micronsPerDegree, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', centerPSF);

        if (isempty(oi))
            error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
        end
end

% Convert spatial units to microns
support.psfX = psfSupportMinutesX * (micronsPerDegree/60);
support.psfY = psfSupportMinutesY * (micronsPerDegree/60);
support.wave = psfSupportWavelength; 

% Returns the central position in case it was the default.
support.position = position;

% This is an important assumption, so we carry it along.
support.micronsPerDegree = micronsPerDegree;

end
