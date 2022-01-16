function [oi, psf, support, zCoeffs, subjID] = oiPosition(zCoeffDatabase, varargin)
% Return the oi appropriate for a visual field position
%
% Synopsis
%   [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase, varargin);
%
% Brief description
%   We compute an OI that is appropriate for different visual field
%   positions using this method.
%
% Inputs
%   zCoeffDatabase   = {'Polans2015','Artal2012'}
%
% Optional key/val pairs
%   subjectRank      = 3;
%   whichEye         ='right eye';
%   positionDegs     = [-10 0];
%   pupilDiamMM      = 3.0;%
%   center psf       = true;
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
 zCoeffDatabase   = 'Polans2015';
 subjectRank = 3;
 eyeside          = 'right';
 positionDegs     = [0 4];
 pupilDiamMM      = 3.0;
 centerpsf        = true;

 [oi, psf, support, zCoeffs, subjID]  = ...
    oiPosition(zCoeffDatabase, 'position',positionDegs, ...
      'pupil diameter', pupilDiamMM, 'subject rank', subjectRank, ...
      'eye side', eyeside,'center psf',false);

 ieNewGraphWin;
 idx = find(support.w == 650);
 imagesc(support.x, support.y, squeeze(psf(:,:,idx)));
 axis 'square'; colormap(gray); xlabel('arc min'); ylabel('arc min');
 set(gca, 'FontSize', 14); set(gca,'xlim',[-5 5],'ylim',[-5 5]);
 title(sprintf('subject #%d', subjID));

 scene = sceneCreate('gridlines',256,64);
 scene = sceneSet(scene,'fov',1);
 oi = oiCompute(oi,scene);
 % oi = oiCrop(oi,'border');
 oiWindow(oi);

%}
%{
% Simpler, with default values
 [oi, psf, support, zCoeffs, subjID]  = oiPosition('Artal2012');
 ieNewGraphWin;
 idx = find(support.w == 650);
 imagesc(support.x, support.y, squeeze(psf(:,:,idx)));
 axis 'square'; colormap(gray); xlabel('arc min'); ylabel('arc min');
 set(gca, 'FontSize', 14);
 title(sprintf('subject #%d', subjID));
%}


%% Read parameters

varargin = ieParamFormat(varargin);

p = inputParser;
validData = {'Polans2015','Artal2012'};
p.addRequired('zCoeffDatabase',@(x)(ismember(x,validData)));

p.addParameter('position',[0 0], @isvector);
p.addParameter('pupildiameter',3,@isscalar);
p.addParameter('subjectrank',1,@isscalar);
p.addParameter('wave',(400:10:750),@isvector);

% See below
% p.addParameter('correctrefraction',true,@islogical);

validEye = {'left','right'};
p.addParameter('eyeside','right',@(x)(ismember(x,validEye)));

p.addParameter('centerpsf',true,@islogical);

p.parse(zCoeffDatabase,varargin{:});

positionDegs = p.Results.position;
pupilDiamMM  = p.Results.pupildiameter;
subjectRank  = p.Results.subjectrank;
eyeSide      = p.Results.eyeside;  
eyeSide = [eyeSide,' eye'];
centerPSF        = p.Results.centerpsf;
wave             = p.Results.wave;

% We could over-ride the specification in the future
% subtractCentralRefraction = p.Results.correctrefraction;


% Hard parameters.  Hmm.  We can add these to the parameters
micronsPerDegree = 290;
wavefrontSpatialSamples = 201;

%%
switch (zCoeffDatabase)
    case 'Polans2015'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = PolansOptics.constants.subjectRanking;
        subjID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRank]));
        
        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(subjID);
        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            PolansOptics.oiForSubjectAtEccentricity(subjID, ...
            eyeSide, positionDegs, pupilDiamMM, wave, micronsPerDegree, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', centerPSF);
        
        support.x = psfSupportMinutesX;
        support.y = psfSupportMinutesY;
        support.w = psfSupportWavelength;
        
        
    case 'Artal2012'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSujectIDs = ArtalOptics.constants.subjectRanking(eyeSide);
        subjID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRank]));
        
        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(eyeSide, subjID);
        if (positionDegs(2) ~= 0)
            fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
            positionDegs(2) = 0;
        end
        
        [oi, psf, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
            ArtalOptics.oiForSubjectAtEccentricity(subjID, ...
            eyeSide, positionDegs, pupilDiamMM, wave, micronsPerDegree, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', centerPSF);
        
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
