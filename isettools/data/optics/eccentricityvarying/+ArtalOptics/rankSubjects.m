function rankSubjects(eyeForRankAnalysis, gradePSFs)
%
% ArtalOptics.rankSubject(true)
%
	assert(ismember(eyeForRankAnalysis, {'right eye', 'left eye'}), 'eye must be set to either ''right eye'', or ''left eye''.');

	% Which PSFs to use for ranking.
	% Here, use the PSFs from the right eye at the fovea (i.e., at eccentricity = 0,0)
	eccentricityDegsPositionForRankAnalysis = struct(...
		'x', 0, ...
		'y', 0);

	switch (eyeForRankAnalysis)
		case 'left eye'
			analyzedSubjectIDs = sort([ ...
                          34    91    85    90    69    15 ...
                          33     7    46    60   130   128 ...
                          21    24    29    35    39    16 ...
                           2     8    73    37   102   105 ...
                          20    68     1    26    41    84 ...
                          80    45    95    62    13    82 ...
                          19   108    57    18    27     4 ...
                          52   119    28    48    97     5  ...
                          127    63    31    56    78    32]); 
        case 'right eye'
        	analyzedSubjectIDs = sort([ ...
                         49    69   101    26    51    70  ...
                         56    65     7    15    81   102  ...
                         68     1    11   129   119    12  ...
                         88   110   111    82     8    22  ...
                         85    97   117    43     2    24  ...
                         45    41   125    99    37    80  ...
                         62    27   103   100    74     9  ...
                         18    34    31     6    84    14  ...
                         48    53    61    44    47    78  ...
                         32]);
    end


    for iSubj = 1:numel(analyzedSubjectIDs)
    	subjectRequiresCentralRefractionCorrection(iSubj) = ...
    		ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(eyeForRankAnalysis, analyzedSubjectIDs(iSubj));
   	end

        

	% Assemble struct with all parameters
	analyzedOpticsParams = struct(...
        'zernikeDataBase',  'Artal2012', ...
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