function [oi, psf, support, zCoeffs, subjID]  = oiPosition(zCoeffDatabase,varargin)
%
% 

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