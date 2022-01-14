function [oi, psf, support, zCoeffs, subjID] = oiPosition(zCoeffDatabase, positionDegs, pupilDiamMM, subjectRankOrder, whichEye)
% Return the oi appropriate for a visual field position
%
% Synopsis
%   [oi, psf, support, zCoeffs,subjID]  = oiPosition(zCoeffDatabase, positionDegs,  pupilDiamMM, subjectRankOrder, whichEye);
%
% Brief description
%   We compute an OI that is appropriate for different visual field
%   positions using this method.
%
% Inputs
%   zCoeffDatabase   = {'Polans2015','Artal2012'}
%   subjectRankOrder = 3;
%   whichEye         ='right eye';
%   positionDegs     = [-10 0];
%   pupilDiamMM      = 3.0;
%
% Optional key/val pairs
%    N/A
%
% Returns
%   oi - an Optical image that can be used to compute for a local scene
%   psf - the PSF of the oi
%   support - Spatial support (minutes of arc)
%                support.x and support.y
%             Wavelength support (nm)
%                support.w
%   zCoeffs - Zernike polynomial coefficients for the oi
%
% See also
%   Was eccDependentOpticsDemo in NC email
%

% Examples:
%{
% Plot the PSF
 zCoeffDatabase   = 'Polans2015';
 subjectRankOrder = 3;
 whichEye         = 'right eye';
 positionDegs     = [0 10];
 pupilDiamMM      = 3.0;

 [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase, positionDegs , pupilDiamMM, subjectRankOrder, whichEye);
 ieNewGraphWin;
 idx = find(support.w == 550);
 imagesc(support.x, support.y, squeeze(psf(:,:,idx)));
 axis 'square'; colormap(gray); xlabel('arc min'); ylabel('arc min');
 set(gca, 'FontSize', 14);
 title(sprintf('%s, subject #%d', zCoeffDatabase, subjID));
%}
%{
 [oi, psf, support, zCoeffs, subjID]  = oiPosition();
 ieNewGraphWin;
 idx = find(support.w == 550);
 imagesc(support.x, support.y, squeeze(psf(:,:,idx)));
 axis 'square'; colormap(gray); xlabel('arc min'); ylabel('arc min');
 set(gca, 'FontSize', 14);
 title(sprintf('%s, subject #%d', zCoeffDatabase, subjID));
%}


%% Read parameters

if notDefined('zCoeffDatabase'), zCoeffDatabase = 'Polans2015'; end
if notDefined('positionDegs'), positionDegs = [0 0]; end
if notDefined('pupilDiamMM'), pupilDiamMM = 3.0; end
if notDefined('subjectRankOrder'), subjectRankOrder = 1; end
if notDefined('whichEye'), whichEye = 'right eye'; end


% Hard parameters.  Hmm.  I guess we need to call this from an object that
% has this information?

wave = 400:10:750;
micronsPerDegree = 290;
zeroCenterPSF = true;
wavefrontSpatialSamples = 201;

%%
switch (zCoeffDatabase)
    case 'Polans2015'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = PolansOptics.constants.subjectRanking;
        subjID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRankOrder]));
        
        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(subjID);
        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            PolansOptics.oiForSubjectAtEccentricity(subjID, ...
            whichEye, positionDegs, pupilDiamMM, wave, micronsPerDegree, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', zeroCenterPSF);
        
        support.x = psfSupportMinutesX;
        support.y = psfSupportMinutesY;
        support.w = psfSupportWavelength;
        
        
    case 'Artal2012'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
        subjID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRankOrder]));
        
        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, subjID);
        if (positionDegs(2) ~= 0)
            fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
            positionDegs(2) = 0;
        end
        
        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            ArtalOptics.oiForSubjectAtEccentricity(subjID, ...
            whichEye, positionDegs, pupilDiamMM, wave, micronsPerDegree, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', zeroCenterPSF);
        
        if (isempty(oi))
            error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
        end
        
        support.x = psfSupportMinutesX;
        support.y = psfSupportMinutesY;
        support.w = psfSupportWavelength;
        
    otherwise
        error('Unknown optics database: ''%s''.', zCoeffDatabase);
end

end
