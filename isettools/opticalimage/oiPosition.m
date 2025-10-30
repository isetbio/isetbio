function [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase,varargin)
% Return oi and many related parameters 
% 

%Example:
%{

%}
%{
    [oi, psf, support, zCoeffs, subjID]  = ...
        oiPosition(zCoeffDatabase, 'position',positionDegs, ...
        'pupil diameter', pupilDiamMM, 'subject rank', thisSubject, ...
        'wave',wave, ...
        'eye side', eyeside,'center psf',centerpsf);
%}

%% Parse parameters

varargin = ieParamFormat(varargin);
p = inputParser;
p.addParameter('zCoeffDatabase', @(x)(ismember(x, {'Polans2015', 'Artal2012'})));
p.addParameter('position',[0 3],@isvector);
p.addParameter('pupildiameter',3,@isnumeric)
p.addParameter('subjectrank',1,@isnumeric);
p.addParameter('wave',550,@isscalar);
p.addParameter('eyeside','right',@(x)(ismember(x,{'left','right'})));
p.addParameter('centerpsf',true,@islogical);
    
p.parse(zCoeffDatabase,varargin{:});

%%
switch (zCoeffDatabase)
    case 'Polans2015'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSubjectIDs = PolansOptics.constants.subjectRanking;
        testSubjectID = rankedSubjectIDs(subjectRankOrder);
        
        % Determine if we need to subtract the subject's central refraction to
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);
        
    case 'Artal2012'
        % Obtain subject IDs ranking in decreasing foveal resolution
        rankedSubjectIDs = ArtalOptics.constants.subjectRanking(whichEye);
        testSubjectID = rankedSubjectIDs(subjectRankOrder);
        
        % Determine if we need to subtract the subject's central refraction to
        subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
end



end

%{
function eccDependentOptics(opticsZernikeCoefficientsDataBase, subjectRankOrder, whichEye, targetEcc, pupilDiamMM )
% 
%{
  USAGE:
       subjectRankOrder = 3;
       whichEye ='right eye';
       targetEcc = [-10 0];
       pupilDiamMM = 3.0;
       opticsZernikeCoefficientsDataBase = 'Polans2015';
       eccDependentOpticsDemo('Polans2015', subjectRankOrder, whichEye, xyEccDegs , pupilDiamMM);
       eccDependentOpticsDemo('Artal2012',  subjectRankOrder, whichEye, xyEccDegs , pupilDiamMM);
%}

   
    wave = 400:10:750;
    micronsPerDegree = 290;
    zeroCenterPSF = true;
    wavefrontSpatialSamples = 201;

    switch (opticsZernikeCoefficientsDataBase)
        case 'Polans2015'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = PolansOptics.constants.subjectRanking;
            testSubjectID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRankOrder]));
            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);
            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
                PolansOptics.oiForSubjectAtEccentricity(testSubjectID, ...
                    whichEye, xyEccDegs, pupilDiamMM, wave, micronsPerDegree, ...
                    'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'zeroCenterPSF', zeroCenterPSF);
               
        case 'Artal2012'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
            testSubjectID = rankedSujectIDs(min([numel(rankedSujectIDs) subjectRankOrder]));
            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
            if (targetEcc(2) ~= 0)
                fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
                targetEcc(2) = 0;
            end
            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
                ArtalOptics.oiForSubjectAtEccentricity(testSubjectID, ...
                    whichEye, targetEcc, pupilDiamMM, wave, micronsPerDegree, ...
                    'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'zeroCenterPSF', zeroCenterPSF);
               
             if (isempty(theOI))
                    error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
             end
        otherwise
            error('Unknown optics database: ''%s''.', opticsZernikeCoefficientsDataBase);
    end

    % Plot the PSF
    hFig = figure();
    set(hFig, 'Color', [1 1 1]);
    idx = find(psfSupportWavelength == 550);
    imagesc(psfSupportMinutesX, psfSupportMinutesY, squeeze(thePSF(:,:,idx)));
    axis 'square'
    colormap(gray);
    xlabel('arc min');
    ylabel('arc min');
    set(gca, 'FontSize', 14);
    title(sprintf('%s, subject #%d', opticsZernikeCoefficientsDataBase, testSubjectID));
   
end
%}