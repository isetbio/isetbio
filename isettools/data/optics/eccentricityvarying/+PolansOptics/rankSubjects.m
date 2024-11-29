function rankSubjects(gradePSFs)
%
% PolansOptics.rankSubject(true)
%

	% Which PSFs to use for ranking.
	% Here, use the PSFs from the right eye at the fovea (i.e., at eccentricity = 0,0)
	eyeForRankAnalysis = 'right eye';
	eccentricityDegsPositionForRankAnalysis = struct(...
		'x', 0, ...
		'y', 0);

    analyzedSubjectIDs = 1:10;
    for iSubj = 1:numel(analyzedSubjectIDs)
    	subjectRequiresCentralRefractionCorrection(iSubj) = ...
    		PolansOptics.constants.subjectRequiresCentralRefractionCorrection(analyzedSubjectIDs(iSubj));
    end

	% Assemble struct with all parameters
	analyzedOpticsParams = struct(...
        'zernikeDataBase',  'Polans2015', ...
        'subjectIndices', analyzedSubjectIDs, ...
        'subjectRequiresCentralRefractionCorrection', subjectRequiresCentralRefractionCorrection, ...
        'zeroCenterPSF', true, ...
        'flipPSFUpsideDown', true, ...
        'pupilDiameterMM', 3, ...
        'whichEye', eyeForRankAnalysis, ...
        'eccentricityPosDegs', eccentricityDegsPositionForRankAnalysis);

	% Where the analyzed data live
	exportsDir = fullfile(isetbioRootPath(), 'calculators/opticsAssessment/exports');

	% Go !
	if (gradePSFs)
		rankingOptics.gradePSFs(exportsDir, analyzedOpticsParams)
		rankedSubjectIDs = rankingOptics.rankGradedPSFs(exportsDir, analyzedOpticsParams);
	else
		rankedSubjectIDs = rankingOptics.rankGradedPSFs(exportsDir, analyzedOpticsParams);
	end
end