function rankSubjects(eyeForRankAnalysis, gradePSFs)
%
% ThibosOptics.rankSubject(true)
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
                          27    66    67     2    26    55     1  ...
                          62    49    43    34    25    28    24  ...
                          17    39    21    52    10    16    63  ...
                          33    38    59    22    42    53     5  ...
                          11    14     4    35    19    56     6  ...
                          37    40    60    29    30    18    44  ...
                          32    20    31    36    48    65     8  ...
                          50    15    46    45    58    41     7  ...
                          13    61    47     9    54    23    64  ...
                           3    57    12    51    69    68    70]); 
        case 'right eye'
        	analyzedSubjectIDs = sort([ ...
                         26    62    43    66     5    31    47  ...
                          14    38    34    48    21    40    50  ...
                           3    58    35    17    16    59    27  ...
                          65    42    41    63    49     8    36  ...
                           1    24    33    55    29    22     2  ...
                          10    28    11    53    18    19    57  ...
                          52    32    25     4     7    54     6  ...
                          60    45    67    37    12    39    56  ...
                          30    51    15    46    20    44     9  ...
                          23    13    61    64    68    69    70]);
    end


    for iSubj = 1:numel(analyzedSubjectIDs)
    	subjectRequiresCentralRefractionCorrection(iSubj) = ...
    		ThibosOptics.constants.subjectRequiresCentralRefractionCorrection(eyeForRankAnalysis, analyzedSubjectIDs(iSubj));
   	end

        

	% Assemble struct with all parameters
	analyzedOpticsParams = struct(...
        'zernikeDataBase',  'Thibos2002', ...
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